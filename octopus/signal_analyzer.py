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
from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
from .npinterface import get_calibration
from weakref import proxy
from itertools import groupby, takewhile
import numpy as np
import json
from collections import deque
from scipy.signal import medfilt

__all__ = ['SignalAnalyzer', 'SignalAnalysis']


# Index numbers in fast5 event tables
EVENT_START_INDEX = 1
EVENT_MOVE_INDEX = 5
EVENT_KMER_SIZE = 5


class SignalAnalyzer:

    def __init__(self, config):
        self.config = config
        self.partmodel = load_partitioner_model(config['partitioner_model'])
        self.sigdump_file = config['sigdump_file']

    def process(self, filename):
        siganalysis = SignalAnalysis(filename, self)
        siganalysis.process()
        return siganalysis


class SignalAnalysis:

    def __init__(self, filename, analyzer):
        self.filename = filename
        self.config_headproc = analyzer.config['head_signal_processing']
        self.analyzer = proxy(analyzer)
        self.error_flags = set()
        self.f5 = None
        self.open_fast5(filename)

    def __del__(self):
        if self.f5 is not None:
            self.f5.close()
            self.f5 = None

    def open_fast5(self, filename):
        self.f5 = Fast5File(filename, 'r')
        self.sampling_rate = self.f5.get_channel_info()['sampling_rate']

    def load_head_signal(self):
        assert self.f5 is not None

        peek_start = self.config_headproc['peeking_start']
        peek_end = self.config_headproc['peeking_end']
        peek_start = int(self.sampling_rate * peek_start)
        peek_end = int(self.sampling_rate * peek_end)

        # Load raw signal and cut out the peek region.
        sig = self.f5.get_raw_data(scale=True)
        sig = np.array(sig[peek_start:peek_end])

        # Calibrate the signal.
        event_size, scale, shift, drift, var, scale_sd, var_sd = (
            get_calibration(self.filename))
        if event_size == 0: # calibration failed in nanopolish
            self.error_flags.add('not_calibrated')
        else:
            sig = sig - shift
            if not np.isclose(drift, 0.):
                sig -= (np.arange(peek_start, peek_start + len(sig)) /
                        self.sampling_rate) * drift
            sig /= scale

        # Filter noises
        outlierpos = np.where((sig > 200.) | (sig < 40.))[0]
        if len(outlierpos) > 0:
            sig[outlierpos] = 100.
        filter_size = self.config_headproc['median_filter_size']
        sig_filtered = np.array(medfilt(sig, filter_size))
        first_index = peek_start + filter_size // 2

        return sig_filtered, first_index

    def detect_partitions(self):
        assert self.f5 is not None

        try:
            with Basecall1DTools(self.f5) as bcall:
                sequence = bcall.get_called_sequence('template')[1][::-1]
                events = bcall.get_event_data('template')
        except KeyError:
            self.error_flags.add('not_basecalled')
            return

        sig, first_index = self.load_head_signal()

        # Run Viterbi fitting to signal model
        plogsum, statecalls = self.analyzer.partmodel.viterbi(sig)

        # Summarize state transitions
        sigparts = {}
        for state, positions in groupby(enumerate(statecalls[1:]), lambda st: st[1][1].name):
            first = last = next(positions)
            for last in positions:
                pass
            sigparts[state] = (first[0] + first_index, last[0] + first_index + 1)

        # TODO: do some quality controls for the state calls by referring to
        # for the time duration, signal emission issue #2

        if 'rna-adapter' not in sigparts:
            self.error_flags.add('adapter_not_detected')
            return # XXX: there can be more chances to detect barcodes by sequences.

        barcode_start, barcode_end = sigparts['rna-adapter']
        # `events' is a 1D array of 1D array rather than a 2D array. Can't use
        # simple multidimensional subscriptions.
        barcode_seqpos = EVENT_KMER_SIZE // 2 + sum(
            ev[EVENT_MOVE_INDEX]
            for ev in takewhile(lambda ev: ev[EVENT_START_INDEX] < barcode_start,
                                events))

        if not self.error_flags:
            print('>{} seqpos={} bpos={} calls={}'.format(self.filename,
                    barcode_seqpos, barcode_start, sigparts))
            print(sequence[barcode_seqpos:barcode_seqpos+25][::-1])
            if self.analyzer.sigdump_file is not None:
                self.dump_signal_partitions(sig, statecalls[1:])

        #print(''.join(self.error_flags), first_index, sig[:10])

    def process(self):
        self.detect_partitions()

    def dump_signal_partitions(self, sig, statecalls):
        outfile = self.analyzer.sigdump_file

        print(0, 'START', statecalls[0][1].name, sep='\t', file=outfile)

        for level, curstate, nextstate in zip(sig, statecalls, statecalls[1:]):
            print(level, curstate[1].name, nextstate[1].name, sep='\t', file=outfile)


# Internal serialization implementation pomegranate to json does not accurately
# recover the original. Use a custom format here.
def load_partitioner_model(modeldata):
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

