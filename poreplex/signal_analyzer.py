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

from pomegranate import (
    HiddenMarkovModel, GeneralMixtureModel, State, NormalDistribution)
from weakref import proxy
from itertools import groupby
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
import multiprocessing as mp
import numpy as np
import pandas as pd
from hashlib import sha1
import h5py
import os
from scipy.signal import medfilt
from .utils import union_intervals

__all__ = ['SignalAnalyzer', 'SignalAnalysis']


class PipelineHandledError(Exception):
    pass


def load_persistence_store(config, analyzer, storage_name='__poreplex_persistence'):
    import importlib.util as imputil, sys

    if storage_name not in sys.modules:
        storage = {
            'segmodel': load_segmentation_model(config['segmentation_model']),
            'unsplitmodel': load_segmentation_model(config['unsplit_read_detection_model']),
            'kmermodel': pd.read_table(config['kmer_model'], header=0, index_col=0),
        }
        storage['kmersize'] = len(storage['kmermodel'].index[0])

        if config['barcoding']:
            from .barcoding import BarcodeDemultiplexer
            storage['demuxer'] = BarcodeDemultiplexer(config['demultiplexing'])

        if config['albacore_onthefly']:
            from .basecall_albacore import AlbacoreBroker
            storage['albacore'] = AlbacoreBroker(config['albacore_configuration'],
                                                 storage['kmersize'])

        fakespec = imputil.spec_from_file_location('fake', 'fake.py')
        persmod = imputil.module_from_spec(fakespec)
        sys.modules[storage_name] = persmod
        setattr(persmod, 'storage', storage)
    else:
        storage = sys.modules[storage_name].storage

    analyzer.segmodel = storage['segmodel']
    analyzer.unsplitmodel = storage['unsplitmodel']
    analyzer.kmermodel = storage['kmermodel']
    analyzer.kmersize = storage['kmersize']
    if 'demuxer' in storage:
        analyzer.demuxer = storage['demuxer']
    if 'albacore' in storage:
        analyzer.albacore = storage['albacore']


class SignalAnalyzer:

    _EVENT_DUMP_FIELD_NAMES = [
        'mean', 'start', 'stdv', 'length', 'model_state',
        'move', 'p_model_state', 'weights', 'pos', 'end', 'scaled_mean']
    _EVENT_DUMP_FIELD_DTYPES = [
        '<f4', '<u8', '<f4', '<u8', None, '<i4',
        '<f4', '<f4', '<u8', '<u8', '<f8']
    EVENT_DUMP_FIELDS = list(zip(_EVENT_DUMP_FIELD_NAMES, _EVENT_DUMP_FIELD_DTYPES))

    def __init__(self, config, batchid):
        load_persistence_store(config, self)

        self.config = config
        self.inputdir = config['inputdir']
        self.outputdir = config['outputdir']
        self.workerid = sha1(mp.current_process().name.encode()).hexdigest()[:16]
        self.batchid = batchid
        self.formatted_batchid = format(batchid, '08d')
        self.EVENT_DUMP_FIELDS[4] = (self.EVENT_DUMP_FIELDS[4][0], 'S{}'.format(self.kmersize))

        if config['dump_adapter_signals']:
            self.adapter_dump_file, self.adapter_dump_group = \
                self.open_dump_file('adapter-dumps', 'adapter')
            self.adapter_dump_list = []
        else:
            self.adapter_dump_file = self.adapter_dump_group = None

        if config['dump_basecalls']:
            self.basecall_dump_file, self.basecall_dump_group = \
                self.open_dump_file('events', 'basecalled_events')
        else:
            self.basecall_dump_file = self.basecall_dump_group = None

    def process(self, filename):
        try:
            return SignalAnalysis(filename, self).process()
        except Exception as exc:
            return {
                'filename': filename,
                'status': 'unknown_error',
                'error_message': '{}: {}'.format(type(exc).__name__, str(exc))
            }

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

    def predict_barcode_labels(self):
        return self.demuxer.predict()

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

    def __init__(self, filename, analyzer):
        self.filename = filename
        self.config = analyzer.config
        self.analyzer = proxy(analyzer)
        self.fast5 = None
        self.sequence = None
        self.open_data_files(filename)

    def __enter__(self):
        return self

    def __del__(self):
        if self.fast5 is not None:
            self.fast5.close()
        return False

    def open_data_files(self, filename):
        fast5path = os.path.join(self.config['inputdir'], filename)
        self.fast5 = Fast5File(fast5path, 'r')
        if len(self.fast5.status.read_info) != 1:
            raise ValueError('Zero or 2+ reads found in a FAST5 file.')

        self.read_info = self.fast5.status.read_info[0]
        self.channel_info = self.fast5.get_channel_info()
        tracking_id = self.fast5.get_tracking_id()

        self.sampling_rate = self.channel_info['sampling_rate']

        self.metainfo = {
            'filename': filename,
            'read_id': self.read_info.read_id,
            'channel': self.channel_info['channel_number'],
            'start_time': round(self.read_info.start_time / self.sampling_rate, 3),
            'run_id': tracking_id['run_id'],
            'sample_id': tracking_id['sample_id'],
            'duration': self.read_info.duration,
            'num_events': self.read_info.event_data_count,
            'sequence_length': 0,
            'mean_qscore': 0.,
        }

    def load_events(self):
        assert self.fast5 is not None

        if self.config['albacore_onthefly']: # Call albacore to get basecalls.
            events = self.load_events_from_albacore()
        else: # Load from Analysis/ objects in the FAST5.
            events = self.load_events_from_fast5()

        events['pos'] = np.cumsum(events['move'])

        # Rescale the signal to fit in the kmer models
        scaling_params = self.compute_scaling_parameters(events)

        duration = np.hstack((np.diff(events['start']), [1])).astype(np.uint64)
        events['end'] = events['start'] + duration

        events['scaled_mean'] = np.poly1d(scaling_params)(events['mean'])

        return events

    def load_events_from_fast5(self):
        # Load events (15-sample chunks in albacore).
        with Basecall1DTools(self.fast5) as bcall:
            events = bcall.get_event_data('template')
            if events is None:
                raise PipelineHandledError('not_basecalled')

            bcall_summary = self.fast5.get_summary_data(bcall.group_name)['basecall_1d_template']
            self.metainfo['sequence_length'] = bcall_summary['sequence_length']
            self.metainfo['mean_qscore'] = bcall_summary['mean_qscore']
            self.metainfo['num_events'] = bcall_summary['num_events']

            self.sequence = bcall.get_called_sequence('template')[1:] + (0,)

            return pd.DataFrame(events)

    def load_events_from_albacore(self):
        rawdata = self.fast5.get_raw_data(scale=True)
        bcall = (
            self.analyzer.albacore.basecall(
                rawdata, self.channel_info, self.read_info,
                os.path.basename(self.filename).rsplit('.', 1)[0]))
        if bcall is None:
            raise PipelineHandledError('not_basecalled')

        self.metainfo['sequence_length'] = bcall['sequence_length']
        self.metainfo['mean_qscore'] = bcall['mean_qscore']
        self.metainfo['num_events'] = bcall['called_events']

        self.sequence = bcall['sequence'], bcall['qstring'], 0
        return bcall['events']

    def compute_scaling_parameters(self, events):
        # Get median value for each kmer state and match with the ONT kmer model.
        events_summarized = events.groupby('pos', sort=False,
                                           group_keys=False, as_index=False).agg(
                                           {'mean': 'median', 'model_state': 'first'})
        if len(events_summarized) < self.config['head_signal_processing']['minimum_kmer_states']:
            raise PipelineHandledError('too_few_events')

        ev_with_mod = pd.merge(events_summarized, self.analyzer.kmermodel,
                               how='left', left_on='model_state', right_index=True)

        # Filter possible outliers out
        meanratio = ev_with_mod['mean'] / ev_with_mod['level_mean']
        inliers = ev_with_mod[meanratio.between(*
            np.percentile(meanratio,
                self.config['head_signal_processing']['scaler_outlier_trim_range']))]

        # Do the final regression
        return np.polyfit(inliers['mean'], inliers['level_mean'], 1)

    def load_raw_signal(self, scaler, start, end):
        end = min(int(self.fast5.status.read_info[0].duration), end)
        raw_sig = self.fast5.get_raw_data(start=start, end=end, scale=True)
        return scaler(medfilt(raw_sig, self.config['head_signal_processing']['median_filter_size']))

    def detect_segments(self, events):
        headsig = events['scaled_mean'][
                    :(events['pos'] <= self.config['head_signal_processing']['segmentation_scan_limit']).sum()]

        # Run Viterbi fitting to signal model
        plogsum, statecalls = self.analyzer.segmodel.viterbi(headsig)

        # Summarize state transitions
        sigparts = {}
        for _, positions in groupby(enumerate(statecalls[1:]),
                                    lambda st: id(st[1][1])):
            first, state = last, _ = next(positions)
            statename = state[1].name
            #statename = statecalls[1 + first].name
            for last, _ in positions:
                pass
            sigparts[statename] = (first, last) # right-inclusive

        return sigparts

    def detect_unsplit_read(self, events, segments):
        # Detect if the read contains two or more adapters in a single read.
        try:
            payload_start = events.iloc[segments['adapter'][1] + 1]['start']
        except IndexError:
            return False # Must be an adapter-only read

        # Bind settings into the local namespace
        config = self.config['head_signal_processing']['unsplit_read_detection']
        _ = lambda name: int(config[name] * self.sampling_rate)
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

    def trim_adapter(self, events, segments):
        if self.sequence is None:
            return

        adapter_end = segments['adapter'][1]
        kmer_lead_size = self.analyzer.kmersize // 2
        adapter_basecall_length = events.iloc[adapter_end]['pos'] + kmer_lead_size

        if adapter_basecall_length > len(self.sequence[0]):
            raise PipelineHandledError('basecall_table_incomplete')
        elif adapter_basecall_length > 0:
            self.sequence = self.sequence[0], self.sequence[1], adapter_basecall_length

    def push_barcode_signal(self, events, segments):
        adapter_events = events.iloc[slice(*segments['adapter'])]
        if len(adapter_events) > 0:
            signal = adapter_events['scaled_mean'].values
            self.analyzer.demuxer.push(self.metainfo['read_id'], signal)

    def process(self):
        error_set = 'okay'

        try:
            events = self.load_events()
            if self.config['dump_basecalls']:
                self.analyzer.write_basecalled_events(
                        self.metainfo['read_id'], events)

            segments = self.detect_segments(events)
            if 'adapter' not in segments:
                raise PipelineHandledError('adapter_not_detected')

            if self.config['trim_adapter']:
                self.trim_adapter(events, segments)

            if self.config['filter_unsplit_reads']:
                isunsplit_read = self.detect_unsplit_read(events, segments)
                if isunsplit_read:
                    raise PipelineHandledError('unsplit_read')

            outname = 'pass'
            if self.config['barcoding']:
                self.push_barcode_signal(events, segments)

        except PipelineHandledError as exc:
            outname = 'artifact' if exc.args[0] in ('unsplit_read',) else 'fail'
            error_set = exc.args[0]
        else:
            if self.config['dump_adapter_signals']:
                self.dump_adapter_signal(events, segments)

        self.metainfo.update({
            'status': error_set,
            'label': outname,
            'fastq': self.sequence,
        })
        return self.metainfo

    def dump_adapter_signal(self, events, segments):
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


# Internal serialization implementation pomegranate to json does not accurately
# recover the original. Use a custom format here.
def load_segmentation_model(modeldata):
    model = HiddenMarkovModel('model')

    states = {}
    for s in modeldata:
        if len(s['emission']) == 1:
            emission = NormalDistribution(*s['emission'][0][:2])
        else:
            weights = np.array([w for _, _, w in s['emission']])
            dists = [NormalDistribution(mu, sigma)
                     for mu, sigma, _ in s['emission']]
            emission = GeneralMixtureModel(dists, weights=weights)
        state = State(emission, name=s['name'])

        states[s['name']] = state
        model.add_state(state)
        if 'start_prob' in s:
            model.add_transition(model.start, state, s['start_prob'])

    for s in modeldata:
        current = states[s['name']]
        for nextstate, prob in s['transition']:
            model.add_transition(current, states[nextstate], prob)

    model.bake()

    return model

