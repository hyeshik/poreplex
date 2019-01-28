#
# Copyright (c) 2018 Institute for Basic Science
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

from weakref import proxy
from itertools import groupby
import multiprocessing as mp
import numpy as np
import pandas as pd
from hashlib import sha1
from io import StringIO
import traceback
import h5py
import sys
import os
from scipy.signal import medfilt
from .worker_persistence import WorkerPersistenceStorage
from .utils import union_intervals

__all__ = ['SignalAnalyzer', 'SignalAnalysis', 'process_batch']


class SignalAnalysisError(Exception):
    pass


# This function must be picklable.
def process_batch(batchid, reads, config):
    try:
        with SignalAnalyzer(config, batchid) as analyzer:
            return analyzer.process(reads)
    except Exception as exc:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        filename = os.path.split(exc_tb.tb_frame.f_code.co_filename)[-1]
        errorf = StringIO()
        traceback.print_exc(file=errorf)

        return (-1, '[{filename}:{lineno}] Unhandled exception {name}: {msg}'.format(
                        filename=filename, lineno=exc_tb.tb_lineno,
                        name=type(exc).__name__, msg=str(exc)), errorf.getvalue())


class SignalAnalyzer:

    _EVENT_DUMP_FIELD_NAMES = [
        'mean', 'start', 'stdv', 'length', 'model_state',
        'move', 'weights', 'pos', 'end', 'scaled_mean']
    _EVENT_DUMP_FIELD_DTYPES = [
        '<f4', '<u8', '<f4', '<u8', None, '<i4',
        '<f4', '<u8', '<u8', '<f8']
    EVENT_DUMP_FIELDS = list(zip(_EVENT_DUMP_FIELD_NAMES, _EVENT_DUMP_FIELD_DTYPES))

    def __init__(self, config, batchid):
        WorkerPersistenceStorage(config).retrieve_objects(self)

        self.config = config
        self.inputdir = config['inputdir']
        self.outputdir = config['outputdir']
        self.workerid = sha1(mp.current_process().name.encode()).hexdigest()[:16]
        self.batchid = batchid
        self.formatted_batchid = format(batchid, '08d')
        self.open_dumps()

    def process(self, reads):
        inputdir = self.config['inputdir']
        results, loaded = [], []

        # Initialize processors and preload signals from fast5
        nextprocs = []
        prepare_loading = self.loader.prepare_loading
        for f5file, read_id in reads:
            if not os.path.exists(os.path.join(inputdir, f5file)):
                results.append({'filename': f5file, 'status': 'disappeared'})
                continue

            try:
                npread = prepare_loading(f5file, read_id)
                if npread.is_stopped():
                    results.append(npread.report())
                else:
                    siganal = SignalAnalysis(npread, self)
                    nextprocs.append(siganal)
                    loaded.append(npread)
            except Exception as exc:
                error = self.pack_unhandled_exception(f5file, read_id, exc, sys.exc_info())
                results.append(error)

        # Determine scaling parameters
        self.loader.fit_scalers()

        # Perform the main analysis procedures
        procs, nextprocs = nextprocs, []
        for siganal in procs:
            try:
                if not siganal.is_stopped():
                    siganal.process()
                    nextprocs.append(siganal)
                else:
                    sys.stdout.flush()
            except Exception as exc:
                f5file = siganal.npread.filename
                read_id = siganal.npread.read_id
                error = self.pack_unhandled_exception(f5file, read_id, exc, sys.exc_info())
                siganal.set_error(error)
            finally:
                siganal.clear_cache()

        # Call barcode identities for demultiplexing
        if self.config['barcoding']:
            self.demuxer.predict()

        # Copy the final results
        for npread in loaded:
            results.append(npread.report())

        return results

    def pack_unhandled_exception(self, f5filename, read_id, exc, excinfo):
        exc_type, exc_obj, exc_tb = excinfo
        srcfilename = os.path.split(exc_tb.tb_frame.f_code.co_filename)[-1]
        errorf = StringIO()
        traceback.print_exc(file=errorf)

        errmsg = ('[{srcfilename}:{lineno}] ({f5filename}#{read_id}) Unhandled '
                  'exception {name}: {msg}\n{exc}'.format(
            srcfilename=srcfilename, lineno=exc_tb.tb_lineno,
            f5filename=f5filename, read_id=read_id, name=type(exc).__name__, msg=str(exc),
            exc=errorf.getvalue()))

        return {
            'filename': f5filename,
            'read_id': read_id,
            'status': 'unknown_error',
            'error_message': errmsg,
        }

    def open_dumps(self):
        self.EVENT_DUMP_FIELDS[4] = (self.EVENT_DUMP_FIELDS[4][0], 'S{}'.format(self.kmersize))

        if self.config['dump_adapter_signals']:
            self.adapter_dump_file, self.adapter_dump_group = \
                self.open_dump_file('adapter-dumps', 'adapter')
            self.adapter_dump_list = []
        else:
            self.adapter_dump_file = self.adapter_dump_group = None

        if self.config['dump_basecalls']:
            self.basecall_dump_file, self.basecall_dump_group = \
                self.open_dump_file('events', 'basecalled_events')
        else:
            self.basecall_dump_file = self.basecall_dump_group = None

    def open_dump_file(self, subdir, parentgroup):
        h5filename = os.path.join(self.outputdir, subdir,
                                  'part-' + self.workerid + '.h5')
        h5 = h5py.File(h5filename, 'a')
        h5group = h5.require_group(parentgroup +
                                   '/' + self.formatted_batchid)
        return h5, h5group

    def push_adapter_signal_catalog(self, read_id, adapter_start, adapter_end):
        self.adapter_dump_list.append((read_id, adapter_start, adapter_end))

    def write_basecalled_events(self, read_id, events, attrs):
        dataset = np.empty(len(events), dtype=self.EVENT_DUMP_FIELDS)
        for name, dtype in self.EVENT_DUMP_FIELDS:
            dataset[name] = events[name]
        try:
            self.basecall_dump_group[read_id] = dataset
            objattrs = self.basecall_dump_group[read_id].attrs
            for attrname, attrvalue in attrs:
                objattrs[attrname] = attrvalue
        except RuntimeError:
            if read_id not in self.basecall_dump_group:
                raise

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if self.adapter_dump_file is not None:
            catgrp = self.adapter_dump_file.require_group('catalog/adapter')
            encodedarray = np.array(self.adapter_dump_list,
                dtype=[('read_id', 'S36'), ('start', 'i8'), ('end', 'i8')])
            catgrp.create_dataset(self.formatted_batchid, shape=encodedarray.shape,
                                  data=encodedarray)
            self.adapter_dump_file.close()

        if self.basecall_dump_file is not None:
            self.basecall_dump_file.close()


class SignalAnalysis:

    def __init__(self, npread, analyzer):
        self.npread = npread
        self.config = analyzer.config
        self.analyzer = proxy(analyzer)

    def set_error(self, error):
        self.npread.set_error(error['status'], error['error_message'])

    def is_stopped(self):
        return self.npread.is_stopped()

    def clear_cache(self):
        self.npread.close()

    def process(self):
        stride = self.config['signal_processing']['rough_signal_stride']

        try:
            # Load the raw signals for segmentation and in-read adapter signals
            signal = self.npread.load_signal(pool=stride)

            # Rough segmentation of signal for extracting adapter signals
            segments = self.detect_segments(signal, stride)
            if 'adapter' not in segments:
                raise SignalAnalysisError('adapter_not_detected')

            # Dump adapter signals on request from the command line
            if self.config['dump_adapter_signals']:
                self.dump_adapter_signal(signal, segments, stride)

            # Queue a barcode identification task with signal
            if self.config['barcoding']:
                self.push_barcode_signal(signal, segments)

            # Measure poly(A) tail signals
            if self.config['measure_polya']:
                if 'polya-tail' in segments:
                    rough_range = segments['polya-tail']
                else:
                    rough_range = segments['adapter'][1] + 1, None
                self.analyzer.polyaanalyzer(self.npread, rough_range, stride)

            # Load basecalled events for further jobs working also in base-space
            events = self.load_events()
            if self.config['dump_basecalls']:
                self.analyzer.write_basecalled_events(
                        self.npread.read_id, events,
                        self.get_dump_attributes(segments, stride))

            # Trim adapter sequences referring to the segmentation and events
            if self.config['trim_adapter']:
                self.trim_adapter(events, segments, stride)

            # Search for the pattern of adapter signals inside the reads
            if self.config['filter_unsplit_reads']:
                isunsplit_read = self.detect_unsplit_read(events, segments, stride)
                if isunsplit_read:
                    raise SignalAnalysisError('unsplit_read')

            # Discard short sequences
            if self.npread.sequence is not None:
                readlength = len(self.npread.sequence[0]) - self.npread.sequence[2]
                if readlength < self.config['minimum_sequence_length']:
                    raise SignalAnalysisError('sequence_too_short')

        except SignalAnalysisError as exc:
            outname = 'artifact' if exc.args[0] in ('unsplit_read',) else 'fail'
            self.npread.set_status(exc.args[0], stop=True)
            self.npread.set_label(outname)
        else:
            self.npread.set_label('pass')

    def get_dump_attributes(self, segments, stride):
        attrlist = []

        if self.npread.scaling_params is not None:
            sp_scale, sp_shift = self.npread.scaling_params
            attrlist.append(('signal_scale', sp_scale))
            attrlist.append(('signal_shift', sp_shift))

        if 'adapter' in segments:
            attrlist.append(('adapter_begin', np.uint32(segments['adapter'][0] * stride)))
            attrlist.append(('adapter_end', np.uint32((segments['adapter'][1] + 1) * stride)))

        if self.npread.polya is not None:
            polya = self.npread.polya
            if 'polya-tail' in segments:
                attrlist.append(('polya_end_debug',
                                 np.uint32((segments['polya-tail'][1] + 1) * stride)))
            attrlist.append(('polya_begin', np.uint32(polya['begin'])))
            attrlist.append(('polya_end', np.uint32(polya['end'])))
            attrlist.append(('spikes', repr(polya['spikes']).encode()))

        return attrlist

    def load_events(self):
        if self.config['albacore_onthefly']: # Call albacore to get basecalls.
            events = self.npread.call_albacore(self.analyzer.albacore)
        else: # Load from Analysis/ objects in the FAST5.
            events = self.npread.load_fast5_events()

        if self.npread.scaling_params is None:
            raise Exception('Signal scaling is not available yet.')

        events['scaled_mean'] = np.poly1d(self.npread.scaling_params)(events['mean'])
        events['pos'] = np.cumsum(events['move'])

        duration = np.hstack((np.diff(events['start']), [1])).astype(np.uint64)
        events['end'] = events['start'] + duration

        return events

    def trim_adapter(self, events, segments, elspan):
        sequence = self.npread.sequence
        if sequence is not None:
            return

        adapter_end = segments['adapter'][1] * elspan
        kmer_lead_size = self.analyzer.kmersize // 2
        adapter_events = events[events['start'] <= adapter_end]
        if len(adapter_events) <= 0:
            return

        adapter_basecall_length = adapter_events['move'].sum() + kmer_lead_size

        if adapter_basecall_length > len(sequence[0]):
            raise SignalAnalysisError('basecall_table_incomplete')
        elif adapter_basecall_length > 0:
            self.npread.set_adapter_trimming_length(adapter_basecall_length)

    def detect_segments(self, signal, elspan):
        scan_limit = self.config['segmentation']['segmentation_scan_limit'] // elspan
        if len(signal) > scan_limit:
            signal = signal[:scan_limit]

        # Run Viterbi fitting to signal model
        plogsum, statecalls = self.analyzer.segmodel.viterbi(signal)

        # Summarize state transitions
        sigparts = {}
        for _, positions in groupby(enumerate(statecalls[1:]),
                                    lambda st: id(st[1][1])):
            first, state = last, _ = next(positions)
            statename = state[1].name
            for last, _ in positions:
                pass
            sigparts[statename] = (first, last) # right-inclusive

        return sigparts

    def detect_unsplit_read(self, events, segments, elspan):
        # Detect if the read contains two or more adapters in a single read.
        try:
            payload_start = (segments['adapter'][1] + 1) * elspan
        except (KeyError, IndexError):
            return False # Must be an adapter-only read

        # Bind settings into the local namespace
        config = self.config['unsplit_read_detection']
        _ = lambda name, rate=self.npread.sampling_rate: int(config[name] * rate)
        window_size = _('window_size'); window_step = _('window_step')
        strict_duration = _('strict_duration')
        duration_cutoffs = [
            (_('loosen_full_length'), _('loosen_dna_length')),
            (_('strict_full_length'), _('strict_dna_length'))]

        excessive_adapters = []

        for left in range(payload_start, events.iloc[-1]['end'], window_step):
            evblock = events[events['start'].between(left, left + window_size)]
            if len(evblock) < 1:
                break

            _, statecalls = self.analyzer.unsplitmodel.viterbi(evblock['scaled_mean'])
            leader_start = None

            # Find two contiguous states leaders -> adapter and compute sum of the durations
            for _, positions in groupby(enumerate(statecalls[1:]), lambda st: id(st[1][1])):
                first, state = last, _ = next(positions)
                statename = state[1].name
                if statename not in ('adapter', 'leader-high', 'leader-low'):
                    leader_start = None
                    continue

                # Find the last index of matching state calls
                for last, _ in positions:
                    pass

                if leader_start is None:
                    leader_start = first

                if statename != 'adapter':
                    continue

                adapter_end = int(evblock.iloc[last]['end'])
                leader_start_in_read = int(evblock.iloc[leader_start]['start'])
                total_duration = adapter_end - leader_start_in_read
                adapter_duration = adapter_end - evblock.iloc[first]['start']
                total_cutoff, adapter_cutoff = duration_cutoffs[
                        (leader_start_in_read - payload_start) <= strict_duration]

                if total_duration >= total_cutoff and adapter_duration >= adapter_cutoff:
                    excessive_adapters.append([leader_start_in_read, 1 + adapter_end])

                leader_start = None

        if not excessive_adapters:
            return False

        adapter_intervals = (
            [[0, payload_start]] + union_intervals(excessive_adapters)
            + [[np.inf, np.inf]])
        basequality_cutoff = config['basecount_quality_limit']
        count_high_quality_reads = lambda tbl: (
            (tbl.groupby('pos').agg({'p_model_state': 'max'})['p_model_state']
                > basequality_cutoff).sum() if len(tbl) >= 0 else 0)
        subread_lengths = [
            count_high_quality_reads(events[events['start'].between(left, right)])
            for (_, left), (right, _) in zip(adapter_intervals[0:], adapter_intervals[1:])]

        subread_hq_length_total = sum(subread_lengths[1:])

        if (subread_hq_length_total > config['subread_basecount_limit'] or
                (subread_hq_length_total + 1) / (subread_lengths[0] + 1)
                    > config['subread_baseratio_limit']):
            return True

        return False

    def push_barcode_signal(self, signal, segments):
        adapter_signal = signal[segments['adapter'][0]:segments['adapter'][1]+1]
        if len(adapter_signal) > 0:
            self.analyzer.demuxer.push(self.npread, adapter_signal)

    def dump_adapter_signal(self, signal, segments, stride):
        adapter_signal = signal[segments['adapter'][0]:segments['adapter'][1]+1]
        if len(adapter_signal) > 0:
            read_id = self.npread.read_id
            try:
                self.analyzer.adapter_dump_group.create_dataset(read_id,
                    shape=(len(adapter_signal),), dtype=np.float32,
                    data=adapter_signal)
            except:
                if read_id in self.analyzer.adapter_dump_group:
                    return
                raise

            start_pos = segments['adapter'][0] * stride
            end_pos = (segments['adapter'][1] + 1) * stride

            self.analyzer.push_adapter_signal_catalog(read_id, start_pos, end_pos)

