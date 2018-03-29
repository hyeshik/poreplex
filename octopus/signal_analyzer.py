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

from pomegranate import HiddenMarkovModel
from ont_fast5_api.fast5_file import Fast5File
from .npinterface import get_calibration
import numpy as np
from collections import deque

__all__ = ['SignalAnalyzer']


def median_filter(sequence, window=5):
    wingsize = window // 2
    q = deque()
    r = []

    for s in sequence:
        q.append(s)
        #print(np.median(q))
        if len(q) > wingsize:
            r.append(np.median(q))
        if len(q) >= window:
            q.popleft()

    for i in range(wingsize):
        r.append(np.median(q))
        q.popleft()

    return r


class SignalAnalyzer:

    def __init__(self, filename, config):
        self.filename = filename
        self.config_headproc = config['head_signal_processing']
        self.f5 = None
        self.open_fast5(filename)
        self.partitioned = False

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
        sig = sig - shift
        sig -= (np.arange(peek_start, peek_start + len(sig)) /
                self.sampling_rate) * drift
        sig /= scale

        # Filter outliers
        outlierpos = np.where((sig > 200.) | (sig < 40.))[0]
        if len(outlierpos) > 0:
            sig[outlierpos] = 100.
        filter_size = self.config_headproc['median_filter_size']
        sig_filtered = np.array(median_filter(sig, filter_size))
        first_index = peek_start + filter_size // 2

        print(sig_filtered)
        print(first_index)

    def detect_partitions(self):
        sig = self.load_head_signal()
        #print(inputfile, *get_calibration(inputfile))

