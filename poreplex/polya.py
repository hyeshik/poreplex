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

from .csupport import detect_events
import numpy as np
import pandas as pd
import os

class PolyASignalAnalyzer:

    CONFIG_SLOTS = [
        'refinement_expansion', 'signal_stdv_max', 'event_detection',
        'signal_z_threshold', 'spike_tolerance', 'spike_max_length',
        'spike_weight',
    ]

    def __init__(self, config):
        for name in self.CONFIG_SLOTS:
            setattr(self, name, config[name])

    def __call__(self, npread, segments, stride):
        if 'polya-tail' not in segments:
            return

        raw_signal = npread.load_signal(pool=None, pad=False)

        rough_begin, rough_end = segments['polya-tail']
        rough_begin = max(0, rough_begin * stride - self.refinement_expansion)
        rough_end = min(len(raw_signal), (rough_end + 1) * stride + self.refinement_expansion)
        polya_signal = raw_signal[rough_begin:rough_end]

        events = pd.DataFrame(detect_events(polya_signal, **self.event_detection))
        longest_event = events.loc[events['length'].idxmax()]
        calc_z_against_longest_event = (
                    lambda ev1, ev2=longest_event: (ev1['mean'] - ev2['mean']) / ev2['stdv'])

        events['z'] = events.apply(calc_z_against_longest_event, axis=1).abs()

        polya_events = self.find_best_polya_interval(events)

        # Quality check for the longest event
        if longest_event['stdv'] < self.signal_stdv_max and len(polya_events) > 0:
            polya_begin = polya_events.iloc[0]['start']
            last_event = polya_events.iloc[-1]
            polya_end = last_event['start'] + last_event['length']

            dwell_time = int(polya_events[polya_events['z'] <=
                                          self.signal_z_threshold]['length'].sum())
            spike_positions = np.where(polya_events['z'] > self.signal_z_threshold)[0]
            spikes = [
                (polya_events.iloc[spk]['length'],) +
                tuple(polya_events.iloc[spk-1:spk+2]['mean'].tolist())
                for spk in spike_positions]

            npread.set_polya_tail({
                'begin': int(polya_begin) + rough_begin,
                'end': int(polya_end) + rough_begin,
                'dwell_time': dwell_time / npread.sampling_rate,
                'spikes': spikes,
            })

    def find_best_polya_interval(self, events):
        nullmtx = np.zeros([len(events) + 1, len(events) + 1], dtype=np.int64)

        matching_scores = nullmtx.copy()
        matching_scores[0, 1:] = [
            v if v > 0 else v * self.spike_weight
            for v in ((events['z'] <= self.signal_z_threshold) * 2 - 1) * events['length']]

        spike_scores = nullmtx.copy()
        spike_scores[0, 1:] = [
            1 if ev.z <= self.signal_z_threshold
              else -int(ev.length > self.spike_max_length)
            for _, ev in events.iterrows()]

        for i in range(1, len(events)+1): # i=start+1
            for j in range(1, len(events)+1): # j=end+1 (inclusive)
                if i > j:
                    continue

                matching_scores[i, j] = matching_scores[i, j - 1] + matching_scores[0, j]

                s = (-1 if spike_scores[i, j - 1] < 0
                     else (self.spike_tolerance
                           if spike_scores[0, j] > 0
                           else (spike_scores[i, j - 1] - 1
                                 if spike_scores[0, j] == 0 else -1)))
                spike_scores[i, j] = s

        spike_scores = (spike_scores >=
                        self.spike_tolerance).astype(np.int64)
        final_scores = matching_scores[1:, 1:] * spike_scores[1:, 1:]
        polya_start, polya_end = np.unravel_index(final_scores.argmax(), (len(events),) * 2)

        if final_scores[polya_start, polya_end] <= 0:
            return [] # no poly(A) found

        return events.iloc[polya_start:polya_end+1]

