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
from progressbar import ProgressBar, NullBar, UnknownLength
from concurrent.futures import (
    ProcessPoolExecutor, CancelledError, ThreadPoolExecutor)
from pysam import BGZFile
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

    return dirname, dirs, files


def process_file(inputfile, outputprefix, analyzer):
    result = analyzer.process(inputfile, outputprefix)
    return result


def process_batch(batchid, reads, config):
    with SignalAnalyzer(config, batchid) as analyzer:
        for f5file in reads:
            # =-= XXX temporary
            #if os.path.exists(outputprefix + '.npy'):
            #    continue
            # XXX
            process_file(f5file, 'x', analyzer)

#    def show_memory_usage():
#        usages = open('/proc/self/statm').read().split()
#        print('{:05d} MEMORY total={} RSS={} shared={} data={}'.format(
#                batchid, usages[0], usages[1], usages[2], usages[4]))
#    show_memory_usage()

        return len(reads), analyzer.sequences


class FASTQWriter:

    def __init__(self, outputdir, barcodes):
        self.outputdir = outputdir
        self.barcodes = barcodes

        self.open_streams()

    def open_streams(self):
        self.streams = {
            barcode: BGZFile(self.get_output_path(barcode), 'w')
            for barcode in self.barcodes}

    def close(self):
        for stream in self.streams.values():
            stream.close()

    def get_output_path(self, name):
        return os.path.join(self.outputdir, 'fastq', name + '.fastq.gz')

    def write_sequences(self, seqpacks):
        for barcode, entries in seqpacks.items():
            formatted = ''.join('@{}\n{}\n+\n{}\n'.format(*fields) for fields in entries)
            self.streams[barcode].write(formatted.encode('ascii'))


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

        self.barcodes = ['Undetermined'] # XXX: change this

    def __enter__(self):
        self.fastq_writer = FASTQWriter(self.args.output, self.barcodes)

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

        self.fastq_writer.close()

    def errx(self, message):
        if self.running:
            errprint(message)
            self.stop('ERROR')

    def stop(self, signalname):
        if self.running:
            errprint("\nTermination request in process. Please wait for a moment.")
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
            nprocessed, seqs = await self.run_in_executor_compute(
                process_batch, batchid, files, self.config)
        except CancelledError:
            return

        self.reads_processed += nprocessed
        self.reads_queued -= nprocessed

        try:
            await self.run_in_executor_io(self.fastq_writer.write_sequences, seqs)
        except CancelledError:
            return

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

    async def scan_dir_recursive(self, dirname, top=False):
        if not self.running:
            return

        try:
            errormsg = None
            path, dirs, files = await self.run_in_executor_io(
                scan_dir_recursive_worker, dirname)
        except CancelledError as exc:
            if top: return
            else: raise exc
        except Exception as exc:
            errormsg = str(exc)

        if errormsg is not None:
            return self.errx('ERROR: ' + str(errormsg))

        for filename in files:
            fullpath = os.path.join(path, filename)
            self.queue_processing(fullpath)

        try:
            for subdir in dirs:
                subdirpath = os.path.join(path, subdir)
                await self.scan_dir_recursive(subdirpath)
        except CancelledError as exc:
            if top: return
            else: raise exc

        if top:
            self.flush_jobstack()
            self.scan_finished = True

    async def monitor_progresses(self):
        pbarclass = NullBar if self.config['quiet'] else ProgressBar
        pbar = pbarclass(max_value=UnknownLength, initial_value=0)
        pbar.start()
        pbar_growing = True

        while self.running:
            if pbar_growing and self.scan_finished:
                pbar_growing = False
                pbar = pbarclass(max_value=self.reads_found,
                                 initial_value=self.reads_processed)
                pbar.start()
            else:
                pbar.update(self.reads_processed)

            await asyncio.sleep(0.3)

            if self.scan_finished and self.reads_queued <= 0:
                break

    def terminate_executors(self):
        force_terminate_executor(self.executor_compute)


    @classmethod
    def run(kls, config, args):
        with kls(config, args) as sess:
            monitor_task = sess.loop.create_task(sess.monitor_progresses())

            scanjob = sess.scan_dir_recursive(args.input, top=True)
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

