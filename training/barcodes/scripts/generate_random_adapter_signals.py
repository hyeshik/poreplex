#!/usr/bin/env python3
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

import h5py
from numpy.lib.format import open_memmap
from random import sample, choice, normalvariate, seed
from bisect import bisect
import numpy as np
import pandas as pd
import sys
import os

START_MARK = -10000.
END_MARK = 10000


def load_adapter_list(adapter_source):
    with h5py.File(adapter_source, 'r') as adpfile:
        adapters = pd.DataFrame(adpfile['catalog/adapter'][:])
        adapters['read_id'] = adapters['read_id'].apply(bytes.decode)
        return adapters

def make_target_distribution(adapters, length_histogram_range,
                             length_histogram_bins, target_population, meansig_stride):
    lengths = (adapters['end'] - adapters['start']) // meansig_stride

    counts, edges = np.histogram(lengths, bins=length_histogram_bins,
                                 range=length_histogram_range)
    counts = pd.Series(counts, index=edges[:-1].round().astype(np.int64))

    assg_real = counts / counts.sum() * target_population
    assg = assg_real.astype(np.int64)
    assg_total = assg.sum()
    if assg_total < target_population:
        compensation = (assg_real - assg).nlargest(target_population - assg_total).clip(1, 1)
        compensation = compensation.reindex(assg.index).fillna(0).astype(np.int64)
        assg += compensation

    return assg.tolist()


def load_single_signal_as_fragments(adpfile, read_id, breakdensity, length_min, length_max):
    adsig = adpfile['adapter/{}/{}'.format(read_id[:3], read_id)]
    adsig = np.hstack([[START_MARK], adsig[:], [END_MARK]])

    breakpoints = sample(range(1, len(adsig)), int((len(adsig) - 1) * breakdensity))
    breakpoints = [0] + sorted(breakpoints) + [len(adsig)]

    slices = np.fromiter(zip(breakpoints, breakpoints[1:]),
                         dtype=[('start', np.int32), ('stop', np.int32)])
    slice_sizes = slices['stop'] - slices['start']

    valid_slices = ((slice_sizes >= length_min) & (slice_sizes <= length_max))
    valid_slices = slices[valid_slices]

    for slicerange in valid_slices:
        sigslice = adsig[slicerange['start']: slicerange['stop']]
        yield (sigslice[0], sigslice)


def load_fragment_source(adapter_source, adapters, fragment_random_source_size,
                         fragment_oversampling_rate, fragment_break_density,
                         fragment_length_min, fragment_length_max):

    selected_read_ids = adapters.sample(fragment_random_source_size)['read_id']

    slice_signals = []
    start_slices = []
    internal_slices = []

    with h5py.File(adapter_source, 'r') as adpfile:
        for i, read_id in enumerate(selected_read_ids.tolist() * fragment_oversampling_rate):
            frags = load_single_signal_as_fragments(adpfile, read_id, fragment_break_density,
                        fragment_length_min, fragment_length_max)

            for mast, signal in frags:
                signal_idx = len(slice_signals)
                slice_signals.append(signal)

                if np.isclose(mast, START_MARK):
                    start_slices.append(signal_idx)
                else:
                    internal_slices.append((mast, signal_idx))

            if i % 500 == 0:
                print(i, end=' ')
                sys.stdout.flush()

        print()

    internal_slices.sort()

    return slice_signals, start_slices, internal_slices

def generate_random_signal(slice_signals, start_slices, internal_slices, stitching_noise):
    fragments = []

    start_idx = choice(start_slices)
    fragments.append(slice_signals[start_idx])

    try:
        while not np.isclose(fragments[-1][-1], END_MARK):
            stkey = fragments[-1][-1] + normalvariate(*stitching_noise)
            found = bisect(internal_slices, (stkey,))
            matching_next = internal_slices[found]
            fragments.append(slice_signals[matching_next[1]][1:])

        return np.hstack(fragments)[1:-1]
    except IndexError:
        return []

if __name__ == '__main__':
    import sys
    import yaml

    seed(922)

    config = yaml.load(open('config.yaml'))['decoy_signal_generation']

    meansig_stride = config['mean_signal_stride']
    fragment_oversampling_rate = config['fragment_oversampling_rate']
    fragment_random_source_size = config['fragment_random_source_size']
    fragment_break_density = config['fragment_break_density']
    fragment_length_min = config['fragment_length_min']
    fragment_length_max = config['fragment_length_max']

    stitching_noise = (0, config['stitching_noise_stdev']) # (mu, sigma)

    length_histogram_range = (config['length_histogram_min'], config['length_histogram_max'])
    length_histogram_bins = config['length_histogram_bins']

    target_population_per_round = config['target_population_per_round']
    num_rounds = config['num_rounds']

    output_length = config['output_length']
    output_file = sys.argv[2]

    adapter_source = sys.argv[1]

    print('==> Loading the adapter list')
    adapters = load_adapter_list(adapter_source)

    print('==> Building target distribution of signal lengths')
    target_distribution = make_target_distribution(adapters, length_histogram_range,
                                         length_histogram_bins,
                                         target_population_per_round,
                                         meansig_stride)

    if os.path.exists(output_file):
        os.unlink(output_file)
    out = open_memmap(output_file, 'w+', dtype=np.float32,
                      shape=(num_rounds * target_population_per_round, output_length))

    for round_no in range(num_rounds):
        print()
        print("## Round {}/{} ##".format(round_no + 1, num_rounds))
        print()

        print('==> Loading random fragments from adapter signals...')
        slice_signals, start_slices, internal_slices = \
            load_fragment_source(adapter_source, adapters, fragment_random_source_size,
                                 fragment_oversampling_rate, fragment_break_density,
                                 fragment_length_min, fragment_length_max)

        print('==> Populating random stitches of the fragments...')
        histtogo = list(target_distribution)
        produced_sigs = []

        while len(produced_sigs) < target_population_per_round:
            rsig = generate_random_signal(slice_signals, start_slices, internal_slices,
                                          stitching_noise)
            binno = int((len(rsig) - length_histogram_range[0]) /
                        (length_histogram_range[1] - length_histogram_range[0]) *
                        length_histogram_bins)
            if binno < 0 or binno >= length_histogram_bins or histtogo[binno] <= 0:
                continue

            histtogo[binno] -= 1

            if len(rsig) > output_length:
                rsig = rsig[-output_length:]
            elif len(rsig) < output_length:
                rsig = np.pad(rsig, [output_length - len(rsig), 0], 'constant')

            produced_sigs.append(rsig)
            if len(produced_sigs) % 100 == 0:
                print(len(produced_sigs), end=' ')
                sys.stdout.flush()

        print()

        print('==> Writing signals out...')
        record_no = round_no * target_population_per_round
        out[record_no:record_no + target_population_per_round] = np.array(produced_sigs)
        out.flush()

