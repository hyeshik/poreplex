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
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
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
        return SignalAnalyzer(config, batchid).process(reads)
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
        'move', 'p_model_state', 'weights', 'pos', 'end', 'scaled_mean']
    _EVENT_DUMP_FIELD_DTYPES = [
        '<f4', '<u8', '<f4', '<u8', None, '<i4',
        '<f4', '<f4', '<u8', '<u8', '<f8']
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
        results, contexts = [], {}

        # Initialize processors and preload signals from fast5
        nextprocs = []
        prepare_loading = self.loader.prepare_loading
        for f5file in reads:
            if not os.path.exists(os.path.join(inputdir, f5file)):
                results.append({'filename': f5file, 'status': 'disappeared'})
                continue

            try:
                ctx = prepare_loading(f5file)
                if ctx['status']:
                    early_errors.append(ctx)
                else:
                    read_id = ctx['meta']['read_id']
                    siganal = SignalAnalysis(ctx, self)
                    nextprocs.append(siganal)
                    contexts[read_id] = ctx
            except Exception as exc:
                error = self.pack_unhandled_exception(filename, exc, sys.exc_info())
                results.append(error)

        # Determine scaling parameters
        self.loader.fit_scalers(contexts)

        # Perform the main analysis procedures
        procs, nextprocs = nextprocs, []
        for siganal in procs:
            try:
                siganal.process()
            except Exception as exc:
                error = self.pack_unhandled_exception(filename, exc, sys.exc_info())
                siganal.context['meta'].update(error)
            else:
                nextprocs.append(siganal)

        # Call barcode identities for demultiplexing
        if self.config['barcoding']:
            self.demuxer.predict(contexts)

        # Copy the final results
        for ctx in contexts.values():
            results.append(ctx['meta'])

        return results

    def pack_unhandled_exception(self, filename, exc, excinfo):
        exc_type, exc_obj, exc_tb = excinfo
        filename = os.path.split(exc_tb.tb_frame.f_code.co_filename)[-1]
        errorf = StringIO()
        traceback.print_exc(file=errorf)

        errmsg = '[{filename}:{lineno}] Unhandled exception {name}: {msg}\n{exc}'.format(
            filename=filename, lineno=exc_tb.tb_lineno,
            name=type(exc).__name__, msg=str(exc), exc=errorf.getvalue())

        return {
            'filename': filename,
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

    def write_basecalled_events(self, read_id, events):
        dataset = np.empty(len(events), dtype=self.EVENT_DUMP_FIELDS)
        for name, dtype in self.EVENT_DUMP_FIELDS:
            dataset[name] = events[name]
        try:
            self.basecall_dump_group[read_id] = dataset
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

    def __init__(self, context, analyzer):
        self.context = context
        self.config = analyzer.config
        self.analyzer = proxy(analyzer)

    def process(self):
        error_set = 'okay'

        try:
            if 'load_signal' not in self.context:
                raise SignalAnalysisError('signal_not_loaded')

            # Load the raw signals for segmentation and in-read adapter signals
            signal, elspan = self.context['load_signal']()

            # Rough segmentation of signal for extracting adapter signals
            segments = self.detect_segments(signal, elspan)
            if 'adapter' not in segments:
                raise SignalAnalysisError('adapter_not_detected')

            # Queue a barcode identification task with signal
            outname = 'pass'
            if self.config['barcoding']:
                self.push_barcode_signal(signal, segments)

            # Load basecalled events for further jobs working also in base-space
            events = self.load_events()

            # Trim adapter sequences referring to the segmentation and events
            if self.config['trim_adapter']:
                self.trim_adapter(events, segments, elspan)

            # Search for the pattern of adapter signals inside the reads
            if self.config['filter_unsplit_reads']:
                isunsplit_read = self.detect_unsplit_read(events, segments, elspan)
                if isunsplit_read:
                    raise SignalAnalysisError('unsplit_read')

        except SignalAnalysisError as exc:
            outname = 'artifact' if exc.args[0] in ('unsplit_read',) else 'fail'
            error_set = exc.args[0]
        else:
            pass
            #if self.config['dump_adapter_signals']:
            #    self.dump_adapter_signal(events, segments)

        self.context['meta'].update({
            'status': error_set,
            'label': outname,
        })

    def load_events(self):
        if self.config['albacore_onthefly']: # Call albacore to get basecalls.
            events = self.load_events_from_albacore()
        else: # Load from Analysis/ objects in the FAST5.
            events = self.load_events_from_fast5()

        events['pos'] = np.cumsum(events['move'])

        duration = np.hstack((np.diff(events['start']), [1])).astype(np.uint64)
        events['end'] = events['start'] + duration

        events['scaled_mean'] = np.poly1d(self.context['signal_scaling'])(events['mean'])

        return events

    def load_events_from_fast5(self):
        # Load events (15-sample chunks in albacore).
        fast5 = self.context['fast5']
        metainfo = self.context['meta']

        with Basecall1DTools(fast5) as bcall:
            events = bcall.get_event_data('template')
            if events is None:
                raise SignalAnalysisError('not_basecalled')

            bcall_summary = fast5.get_summary_data(bcall.group_name)['basecall_1d_template']
            metainfo['sequence_length'] = bcall_summary['sequence_length']
            metainfo['mean_qscore'] = bcall_summary['mean_qscore']
            metainfo['num_events'] = bcall_summary['num_events']
            metainfo['sequence'] = bcall.get_called_sequence('template')[1:] + (0,)

            return pd.DataFrame(events)

    def trim_adapter(self, events, segments, elspan):
        if 'sequence' not in self.context['meta']:
            return

        sequence = self.context['meta']['sequence']
        adapter_end = segments['adapter'][1] * elspan
        kmer_lead_size = self.analyzer.kmersize // 2
        adapter_events = events[events['start'] <= adapter_end]
        if len(adapter_events):
            adapter_basecall_length = adapter_events['move'].sum() + kmer_lead_size
        else:
            adapter_basecall_length = 0

        if adapter_basecall_length > len(sequence[0]):
            raise SignalAnalysisError('basecall_table_incomplete')
        elif adapter_basecall_length > 0:
            self.context['meta']['sequence'] = (
                sequence[0], sequence[1], adapter_basecall_length)

    def load_events_from_albacore(self):
        metainfo = self.context['meta']

        rawdata = self.context['load_signal'](pool=False)[0]
        bcall = (
            self.analyzer.albacore.basecall(
                rawdata, self.context['channel_info'], self.context['read_info'],
                os.path.basename(metainfo['filename']).rsplit('.', 1)[0]))
        if bcall is None:
            raise SignalAnalysisError('not_basecalled')

        metainfo['sequence_length'] = bcall['sequence_length']
        metainfo['mean_qscore'] = bcall['mean_qscore']
        metainfo['num_events'] = bcall['called_events']
        metainfo['sequence'] = bcall['sequence'], bcall['qstring'], 0

        return bcall['events']

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
        _ = lambda name, rate=self.context['sampling_rate']: int(config[name] * rate)
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
        adapter_signal = signal[:segments['adapter'][1]+1]
        if len(adapter_signal) > 0:
            self.analyzer.demuxer.push(self.context['meta']['read_id'], adapter_signal)

    def __dump_adapter_signal(self, events, segments):
        adapter_events = events.iloc[slice(*segments['adapter'])]
        if len(adapter_events) > 0:
            try:
                self.analyzer.adapter_dump_group.create_dataset(self.metainfo['read_id'],
                    shape=(len(adapter_events),), dtype=np.float32,
                    data=adapter_events['scaled_mean'])
            except:
                if self.metainfo['read_id'] in self.analyzer.adapter_dump_group:
                    return
                raise
            self.analyzer.push_adapter_signal_catalog(self.metainfo['read_id'],
                adapter_events['start'].iloc[0], adapter_events['end'].iloc[-1])

