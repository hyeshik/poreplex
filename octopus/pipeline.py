#
# Copyright (c) 2018 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import asyncio
import signal
import sys
import os
from errno import EXDEV
from progressbar import ProgressBar, NullBar, UnknownLength
from concurrent.futures import (
    ProcessPoolExecutor, CancelledError, ThreadPoolExecutor)
from . import *
from .io import (
    FASTQWriter, SequencingSummaryWriter, create_adapter_dumps_inventory,
    create_events_inventory)
from .signal_analyzer import SignalAnalyzer
from .utils import *


def force_terminate_executor(executor):
    executor._call_queue.empty()

    if executor._processes:
        alive_pids = set(
            pid for pid, proc in executor._processes.items()
            if proc.is_alive())
        for pid in alive_pids:
            os.kill(pid, signal.SIGKILL)


def scan_dir_recursive_worker(dirname, suffix='.fast5'):
    files, dirs = [], []
    for entryname in os.listdir(dirname):
        if entryname.startswith('.'):
            continue

        fullpath = os.path.join(dirname, entryname)
        if os.path.isdir(fullpath):
            dirs.append(entryname)
        elif entryname.lower().endswith(suffix):
            files.append(entryname)

    return dirs, files


def show_memory_usage():
    usages = open('/proc/self/statm').read().split()
    print('{:05d} MEMORY total={} RSS={} shared={} data={}'.format(
            batchid, usages[0], usages[1], usages[2], usages[4]))

def process_file(inputfile, outputprefix, analyzer):
    result = analyzer.process(inputfile, outputprefix)
    return result

def process_batch(batchid, reads, config):
    try:
        with SignalAnalyzer(config, batchid) as analyzer:
            return [process_file(f5file, 'x', analyzer) for f5file in reads]
    except Exception as exc:
        import traceback
        traceback.print_exc()
        return -1, 'Unhandled exception {}: {}'.format(type(exc).__name__, str(exc))

class ProcessingSession:

    def __init__(self, config, args):
        self.running = True
        self.scan_finished = False
        self.reads_queued = self.reads_found = 0
        self.reads_processed = 0
        self.next_batch_id = 0
        self.jobstack = []

        self.config = config
        self.args = args

        self.executor_compute = ProcessPoolExecutor(args.parallel)
        self.executor_io = ThreadPoolExecutor(2)

        self.loop = self.fastq_writer = None
        self.pbar = None

    def __enter__(self):
        self.fastq_writer = FASTQWriter(self.config['outputdir'],
                                        self.config['output_names'])
        self.seqsummary_writer = SequencingSummaryWriter(self.config['outputdir'])

        self.loop = asyncio.get_event_loop()
        self.executor_compute.__enter__()
        self.executor_io.__enter__()

        for signame in 'SIGINT SIGTERM'.split():
            self.loop.add_signal_handler(getattr(signal, signame),
                                         self.stop, signame)
        return self

    def __exit__(self, *args):
        self.executor_io.__exit__(*args)
        self.executor_compute.__exit__(*args)
        self.loop.close()

        self.seqsummary_writer.close()
        self.fastq_writer.close()

    def errx(self, message):
        if self.running:
            errprint(message)
            self.stop('ERROR')

    def stop(self, signalname='unknown'):
        if self.running:
            errprint("\nTermination in process. Please wait for a moment.")
            self.running = False
        for task in asyncio.Task.all_tasks():
            task.cancel()

        self.loop.stop()

    def run_in_executor_compute(self, *args):
        return self.loop.run_in_executor(self.executor_compute, *args)

    def run_in_executor_io(self, *args):
        return self.loop.run_in_executor(self.executor_io, *args)

    async def run_process_batch(self, batchid, files):
        try:
            results = await self.run_in_executor_compute(
                process_batch, batchid, files, self.config)

            if len(results) > 0 and results[0] == -1: # Unhandled exception occurred
                error_message = results[1]
                errprint("\nERROR: " + error_message)
                self.stop()
                return

            await self.run_in_executor_io(self.fastq_writer.write_sequences, results)

            if self.config['fast5_output']:
                await self.run_in_executor_io(self.link_fast5, results)

            await self.run_in_executor_io(self.seqsummary_writer.write_results, results)

        except CancelledError:
            return

        self.reads_processed += len(results)
        self.reads_queued -= len(results)

    def link_fast5(self, results):
        indir = self.config['inputdir']
        outdir = self.config['outputdir']
        symlinkfirst = self.config['fast5_always_symlink']
        blacklist_hardlinks = set()

        for entry in results:
            original_fast5 = os.path.join(indir, entry['filename'])
            link_path = os.path.join(outdir, 'fast5', entry['label'], entry['filename'])

            original_dir = os.path.dirname(original_fast5)
            link_dir = os.path.dirname(link_path)

            if not os.path.isdir(link_dir):
                os.makedirs(link_dir)

            if not symlinkfirst and (original_dir, link_dir) not in blacklist_hardlinks:
                try:
                    os.link(original_fast5, link_path)
                except OSError as exc:
                    if exc.errno != EXDEV:
                        raise
                    blacklist_hardlinks.add((original_dir, link_dir))
                else:
                    continue

            os.symlink(os.path.abspath(original_fast5), link_path)

    def queue_processing(self, filepath):
        self.jobstack.append(filepath)
        self.reads_queued += 1
        self.reads_found += 1
        if len(self.jobstack) >= self.config['batch_chunk_size']:
            self.flush_jobstack()

    def flush_jobstack(self):
        if self.running and self.jobstack:
            batch_id = self.next_batch_id
            self.next_batch_id += 1
            work = self.run_process_batch(batch_id, self.jobstack[:])
            self.loop.create_task(work)
            del self.jobstack[:]

    async def scan_dir_recursive(self, topdir, dirname=''):
        if not self.running:
            return

        is_topdir = (dirname == '')

        try:
            errormsg = None
            dirs, files = await self.run_in_executor_io(
                scan_dir_recursive_worker, os.path.join(topdir, dirname))
        except CancelledError as exc:
            if is_topdir: return
            else: raise exc
        except Exception as exc:
            errormsg = str(exc)

        if errormsg is not None:
            return self.errx('ERROR: ' + str(errormsg))

        for filename in files:
            filepath = os.path.join(dirname, filename)
            self.queue_processing(filepath)

        try:
            for subdir in dirs:
                subdirpath = os.path.join(dirname, subdir)
                await self.scan_dir_recursive(topdir, subdirpath)
        except CancelledError as exc:
            if is_topdir: return
            else: raise exc

        if is_topdir:
            self.flush_jobstack()
            self.scan_finished = True

    async def monitor_progresses(self):
        pbarclass = NullBar if self.config['quiet'] else ProgressBar
        self.pbar = pbarclass(max_value=UnknownLength, initial_value=0)
        self.pbar.start()
        pbar_growing = True

        while self.running:
            if pbar_growing and self.scan_finished:
                pbar_growing = False
                self.pbar = pbarclass(max_value=self.reads_found,
                                 initial_value=self.reads_processed)
                self.pbar.start()
            else:
                self.pbar.update(self.reads_processed)

            await asyncio.sleep(0.3)

            if self.scan_finished and self.reads_queued <= 0:
                break

    def terminate_executors(self):
        force_terminate_executor(self.executor_compute)

    def finalize_results(self):
        if self.config['dump_adapter_signals']:
            print("==> Creating an inventory for adapter signal dumps...")
            adapter_dump_prefix = os.path.join(self.config['outputdir'], 'adapter-dumps')
            create_adapter_dumps_inventory(
                os.path.join(adapter_dump_prefix, 'inventory.h5'),
                os.path.join(adapter_dump_prefix, 'part-*.h5'))

        if self.config['dump_basecalls']:
            print("==> Creating an inventory for basecalled events...")
            events_prefix = os.path.join(self.config['outputdir'], 'events')
            create_events_inventory(
                os.path.join(events_prefix, 'inventory.h5'),
                os.path.join(events_prefix, 'part-*.h5'))

    @classmethod
    def run(kls, config, args):
        with kls(config, args) as sess:
            print("==> Processing FAST5 files...")
            monitor_task = sess.loop.create_task(sess.monitor_progresses())

            scanjob = sess.scan_dir_recursive(config['inputdir'])
            sess.loop.create_task(scanjob)

            try:
                sess.loop.run_until_complete(monitor_task)
            except CancelledError:
                errprint('\nInterrupted')
            except Exception as exc:
                errprint('\nERROR: ' + str(exc))

            sess.terminate_executors()

            for task in asyncio.Task.all_tasks():
                if not (task.done() or task.cancelled()):
                    try:
                        sess.loop.run_until_complete(task)
                    except CancelledError:
                        errprint('\nInterrupted')
                    except Exception as exc:
                        errprint('\nERROR: ' + str(exc))

            if (not config['quiet'] and
                    sess.scan_finished and sess.reads_found == sess.reads_processed):
                if sess.pbar is not None:
                    sess.pbar.finish()

                print('\n')

                sess.finalize_results()

                print('==> Finished.')
