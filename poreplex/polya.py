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
from scipy.signal import medfilt
from functools import partial
import os

class PolyASignalAnalyzer:

    CONFIG_SLOTS = [
        'refinement_expansion', 'event_detection', 'polya_stdv_max', 'polya_stdv_range',
        'spike_tolerance', 'spike_weight', 'openend_expansion', 'recalibrate_shifted_signal',
        'polya_mean_dist', 'polya_mean_z_cutoff', 'polya_mean_trigger_recalibration',
        'maximum_openend_extension', 'median_pre_filter',
    ]

    def __init__(self, config):
        for name in self.CONFIG_SLOTS:
            setattr(self, name, config[name])

        mean_loc, mean_scale = config['polya_mean_dist']
        self.polya_mean_cutoff = (
            mean_loc - mean_scale * config['polya_mean_z_cutoff'],
            mean_loc + mean_scale * config['polya_mean_z_cutoff'])

        self.polya_mean_trigger_recalibration *= config['polya_mean_dist'][1]

    def __call__(self, npread, rough_range, stride, polya_range=None, ext_depth=0):
        raw_signal = npread.load_signal(pool=None, pad=False)

        minimum_expansion_unit = self.openend_expansion // stride
        rough_begin, rough_end = rough_range
        if rough_end is None or rough_end - rough_begin < minimum_expansion_unit:
            rough_end = rough_begin + minimum_expansion_unit

        insp_begin = max(0, rough_begin * stride - self.refinement_expansion)
        insp_end = min(len(raw_signal), (rough_end + 1) * stride + self.refinement_expansion)
        adapter_end = rough_range[0] * stride - insp_begin
        polya_signal = raw_signal[insp_begin:insp_end]
        if self.median_pre_filter > 1:
            polya_signal = medfilt(polya_signal, self.median_pre_filter)

        events = pd.DataFrame(detect_events(polya_signal, **self.event_detection))
        events['end'] = (events['start'] + events['length']).astype(np.int64)
        events['is_polya'] = events['mean'].between(*(polya_range or self.polya_mean_cutoff))

        entryfunc = (self.call_polya if rough_range[1] is not None
                                     else self.try_recalibrate_shifted_signal)
        entryfunc(npread, events, polya_signal, insp_begin, insp_end,
                  (rough_begin, rough_end), adapter_end, len(raw_signal), stride,
                  polya_range, ext_depth)

    def call_polya(self, npread, events, polya_signal, signal_begin, signal_end,
                   base_range, adapter_end, full_length, stride, polya_range=None,
                   ext_depth=0):
        polya_events = self.find_best_polya_interval(events)

        # Retry it with right-extended interval if poly(A) is not terminated within the interval
        if (len(polya_events) > 0 and polya_events.index[-1] == events.index[-1] and
                    signal_end < full_length and ext_depth < self.maximum_openend_extension):
            return self(npread,
                        (base_range[0], base_range[1] + self.openend_expansion // stride),
                        stride, polya_range, ext_depth + 1)

        # Warp to the recalibration mode if the mean signal level is out of expected range.
        def is_polya_signal_shifted():
            mean_polya_level = ((polya_events['mean'] * polya_events['length']).sum() /
                                polya_events['length'].sum())
            return (
                abs(mean_polya_level - self.polya_mean_dist[0]) >
                self.polya_mean_trigger_recalibration)

        if len(polya_events) == 0 or (polya_range is None and is_polya_signal_shifted()):
            return self.try_recalibrate_shifted_signal(npread, events, polya_signal,
                        signal_begin, signal_end, base_range, adapter_end,
                        full_length, stride, None, ext_depth)

        # Quality check for the longest event
        longest_event = polya_events.loc[polya_events['length'].idxmax()]
        longest_event_stdv = self.calc_internal_polya_stdv(polya_signal, longest_event)

        if longest_event_stdv < self.polya_stdv_max and len(polya_events) > 0:
            polya_begin = polya_events.iloc[0]['start']
            last_event = polya_events.iloc[-1]
            polya_end = last_event['start'] + last_event['length']

            dwell_time = int(polya_events[polya_events['is_polya']]['length'].sum())
            spike_positions = np.where(~polya_events['is_polya'])[0]
            spikes = [
                (polya_events.iloc[spk]['length'],) +
                tuple(polya_events.iloc[spk-1:spk+2]['mean'].tolist())
                for spk in spike_positions]

            npread.set_polya_tail({
                'begin': int(polya_begin) + signal_begin,
                'end': int(polya_end) + signal_begin,
                'dwell_time': dwell_time / npread.sampling_rate,
                'spikes': spikes,
            })
        elif polya_range is None:
            self.try_recalibrate_shifted_signal(npread, events, polya_signal,
                        signal_begin, signal_end, base_range, adapter_end,
                        full_length, stride, None, ext_depth)

    def try_recalibrate_shifted_signal(self, npread, events, polya_signal, signal_begin,
                                  signal_end, base_range, adapter_end, full_length, stride,
                                  polya_range=None, ext_depth=0):
        config = self.recalibrate_shifted_signal
        anchorevents = events[(events['start'] <= adapter_end +
                                    config['max_dist_from_adapter']) &
                              (events['end'] > adapter_end) &
                              (events['stdv'] < config['max_stdv'])]
        if len(anchorevents) == 0:
            return

        polya_mean = ((anchorevents['mean'] * anchorevents['length']).sum() /
                      anchorevents['length'].sum())
        polya_range = (
            polya_mean - self.polya_mean_dist[1] * self.polya_mean_z_cutoff,
            polya_mean + self.polya_mean_dist[1] * self.polya_mean_z_cutoff)

        events['is_polya'] = events['mean'].between(*polya_range)
        if events[events['is_polya']]['length'].sum() >= config['min_length']:
            self.call_polya(npread, events, polya_signal, signal_begin, signal_end,
                            base_range, adapter_end, full_length, stride, polya_range,
                            ext_depth)

    def calc_internal_polya_stdv(self, signal, row):
        length = int(row['length'])
        begin = int(row['start'] + length * self.polya_stdv_range[0])
        end = int(row['start'] + length * self.polya_stdv_range[1])
        return signal[begin:end].std() if end - begin > 2 else np.nan

    def find_best_polya_interval(self, events):
        nullmtx = np.zeros([len(events) + 1, len(events) + 1], dtype=np.int64)

        matching_scores = nullmtx.copy()
        matching_scores[0, 1:] = [
            v if v > 0 else v * self.spike_weight
            for v in (events['is_polya'] * 2 - 1) * events['length']]

        spike_scores = nullmtx.copy()
        spike_scores[0, 1:] = np.choose(events['is_polya'],
            np.array([-events['length'], [1] * len(events)]))

        for i in range(1, len(events)+1): # i=start+1
            for j in range(1, len(events)+1): # j=end+1 (inclusive)
                if i > j:
                    continue

                matching_scores[i, j] = matching_scores[i, j - 1] + matching_scores[0, j]

                s = (-1 if spike_scores[i, j - 1] < 0
                     else (self.spike_tolerance if spike_scores[0, j] > 0
                           else spike_scores[i, j - 1] + spike_scores[0, j]))
                spike_scores[i, j] = s

        spike_scores = (spike_scores > 0).astype(np.int64)
        final_scores = matching_scores[1:, 1:] * spike_scores[1:, 1:]
        polya_start, polya_end = np.unravel_index(final_scores.argmax(), (len(events),) * 2)

        if final_scores[polya_start, polya_end] <= 0:
            return [] # no poly(A) found

        return events.iloc[polya_start:polya_end+1]

