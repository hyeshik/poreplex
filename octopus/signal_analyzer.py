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

from pomegranate import (
    HiddenMarkovModel, GeneralMixtureModel, State, NormalDistribution)
from weakref import proxy
from itertools import groupby, takewhile
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
from statsmodels import api as smapi
from statsmodels.formula.api import ols
import numpy as np
import pandas as pd
import os
from scipy.signal import medfilt

__all__ = ['SignalAnalyzer', 'SignalAnalysis']


# Index numbers in fast5 event tables
EVENT_START_INDEX = 1
EVENT_MOVE_INDEX = 5
EVENT_KMER_SIZE = 5


class SignalAnalyzer:

    def __init__(self, config):
        self.config = config
        self.segmodel = load_segmentation_model(config['segmentation_model'])
        self.kmermodel = pd.read_table(config['kmer_model'], header=0, index_col=0)
        self.kmersize = len(self.kmermodel.index[0])
        self.sigdump_file = config['sigdump_file']

    def process(self, filename, outputprefix):
        print(filename)
        siganalysis = SignalAnalysis(filename, outputprefix, self)
        siganalysis.process()
        return siganalysis


class SignalAnalysis:

    def __init__(self, filename, outputprefix, analyzer):
        self.filename = filename
        self.outputprefix = outputprefix
        self.config = analyzer.config['head_signal_processing']
        self.analyzer = proxy(analyzer)
        self.error_flags = set()
        self.fast5 = None
        self.open_data_files(filename)

    def __enter__(self):
        return self

    def __del__(self):
        if self.fast5 is not None:
            self.fast5.close()
        return False

    def open_data_files(self, filename):
        self.fast5 = Fast5File(filename, 'r')
        self.sampling_rate = self.fast5.get_channel_info()['sampling_rate']

    def load_head_events(self):
        assert self.fast5 is not None

        # Load events (15-sample chunks in albacore).
        with Basecall1DTools(self.fast5) as bcall:
            events = bcall.get_event_data('template')
            if events is None:
                self.error_flags.add('not_basecalled')
                return

            events = pd.DataFrame(events)
            events['pos'] = np.cumsum(events['move'])

        # Rescale the signal to fit in the kmer models
        scaling_params = self.compute_scaling_parameters(events)
        if scaling_params is None:
            return

        duration = np.hstack((np.diff(events['start']), [1]))
        events['end'] = (events['start'] + duration).astype(np.uint64)

        return (events[events['pos'] <= self.config['segmentation_scan_limit']],
                np.poly1d(scaling_params))

    def compute_scaling_parameters(self, events):
        # Get median value for each kmer state and match with the ONT kmer model.
        events_summarized = events.groupby('pos', sort=False,
                                           group_keys=False, as_index=False).agg(
                {'mean': 'median', 'start': 'first', 'model_state': 'first'})
        if len(events_summarized) < self.config['minimum_kmer_states']:
            self.error_flags.add('too_few_kmer_states')
            return

        ev_with_mod = pd.merge(events_summarized, self.analyzer.kmermodel,
                               how='left', left_on='model_state', right_index=True)

        # Filter outliers out
        regr = ols('obs ~ model', data={'obs': ev_with_mod['mean'],
                                        'model': ev_with_mod['level_mean']}).fit()
        oltest = regr.outlier_test()
        inliers = ev_with_mod[oltest['bonf(p)'] > self.config['scaler_outlier_adjp']]
        if len(inliers) < self.config['scaler_minimum_inliers']:
            self.error_flags.add('scaling_failed')
            return

        # Do the final regression
        return np.polyfit(inliers['mean'], inliers['level_mean'], 1)

    def load_raw_signal(self, scaler, start, end):
        end = min(int(self.fast5.status.read_info[0].duration), end)
        raw_sig = self.fast5.get_raw_data(start=start, end=end, scale=True)
        return scaler(raw_sig)

    def detect_segments(self):
        headevents = self.load_head_events()
        if headevents is None:
            return
        headevents, scaler = headevents
        headsig = scaler(headevents['mean'])

        # Run Viterbi fitting to signal model
        plogsum, statecalls = self.analyzer.segmodel.viterbi(headsig)

        # Summarize state transitions
        sigparts = {}
        for state, positions in groupby(enumerate(statecalls[1:]), lambda st: st[1][1].name):
            first = last = next(positions)
            for last in positions:
                pass
            sigparts[state] = (first[0], last[0] + 1)

        # TODO: do some quality controls for the state calls by referring to
        # for the time duration, signal emission issue #2

        if 'rna-adapter' not in sigparts or 'dna-adapter' not in sigparts:
            self.error_flags.add('adapter_not_detected')
            return # XXX: there can be more chances than just detecting them by sequences.

        barcode_start, barcode_end = sigparts['rna-adapter']
        dna_start, dna_end = sigparts['dna-adapter']

        if not self.error_flags:
            outdir = os.path.dirname(self.outputprefix)
            if not os.path.isdir(outdir):
                try:
                    os.makedirs(outdir)
                except OSError:
                    pass

            with open(self.outputprefix + '-dna.npy', 'wb') as eventsout:
                dna_raw_start = headevents.iloc[dna_start]['start']
                dna_raw_end = headevents.iloc[dna_end - 1]['end']
                dna_signal = self.load_raw_signal(scaler, dna_raw_start, dna_raw_end)
                np.save(eventsout, dna_signal)

            with open(self.outputprefix + '-rna.npy', 'wb') as eventsout:
                rna_raw_start = headevents.iloc[barcode_start]['start']
                rna_raw_end = headevents.iloc[barcode_end - 1]['end']
                rna_signal = self.load_raw_signal(scaler, rna_raw_start, rna_raw_end)
                np.save(eventsout, rna_signal)

            if self.analyzer.sigdump_file is not None:
                self.dump_signal_segments(sig, statecalls[1:])

        #print(''.join(self.error_flags), first_index, sig[:10])

    def process(self):
        self.detect_segments()

    def dump_signal_segments(self, sig, statecalls):
        outfile = self.analyzer.sigdump_file

        print(0, 'START', statecalls[0][1].name, sep='\t', file=outfile)

        for level, curstate, nextstate in zip(sig, statecalls, statecalls[1:]):
            print(level, curstate[1].name, nextstate[1].name, sep='\t', file=outfile)


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

