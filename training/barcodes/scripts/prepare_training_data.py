#!/usr/bin/env python3
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

import pandas as pd
import numpy as np
import subprocess as sp
import tempfile
import shutil
import h5py
import glob
import sys
import os
from concurrent import futures

OUTPUT_DTYPE = np.float32

class TemporaryDirectory(object):
    def __init__(self, root='.'):
        self.root = root
        self.path = None

    def __enter__(self):
        self.path = tempfile.mkdtemp(dir=self.root)
        return self

    def __exit__(self, type, value, traceback):
        if self.path is not None:
            shutil.rmtree(self.path)

    def __str__(self):
        return self.path

    def all_files(self):
        return sorted(glob.glob(os.path.join(self.path, 'part*')))

    def merge_into_file(self, outfile):
        infiles = self.all_files()
        if infiles:
            sp.check_call(['cat'] + infiles, stdout=outfile)


def process(signal_file, signal_trim_length, output_path, read_ids):
    with h5py.File(signal_file, 'r') as h5, open(output_path, 'wb') as arrayout:
        sigbuf = []
        siggroup = h5['adapter']

        for read_id in read_ids:
            signal = siggroup['{}/{}'.format(read_id[:3], read_id)][:]
            if len(signal) < signal_trim_length:
                signal = np.pad(signal,
                    (signal_trim_length - len(signal), 0), 'constant')
            elif len(signal) > signal_trim_length:
                signal = signal[-signal_trim_length:]

            sigbuf.append(signal.astype(OUTPUT_DTYPE))

        np.array(sigbuf).tofile(arrayout)

    return len(read_ids)


def main(signal_file, signal_trim_length, catalog_input, output_file,
         parallel=8, chunk_size=2000):

    selreads = pd.read_table(catalog_input)

    with futures.ProcessPoolExecutor(parallel) as executor, \
            TemporaryDirectory() as tmpdir:

        jobs = []
        jobbases = np.arange(int(np.ceil(len(selreads) / chunk_size))) * chunk_size
        for jobbase in jobbases:
            job = executor.submit(process, signal_file, signal_trim_length,
                            '{}/part{:012d}'.format(tmpdir, jobbase),
                            selreads['read_id'].iloc[jobbase:jobbase+chunk_size].tolist())
            jobs.append(job)

        done = 0
        for job in jobs:
            done += job.result()
            print('\r{:,} / {:,} files ({:.2f}%)'.format(
                    done, len(selreads), done / len(selreads) * 100), end='')
            sys.stdout.flush()

        print('\nMerging...')
        tmpdir.merge_into_file(open('{}/merged'.format(tmpdir), 'wb'))

        print('\nConverting...')
        elementsize = signal_trim_length
        arr = np.frombuffer(open('{}/merged'.format(tmpdir), 'rb').read(),
                            dtype=OUTPUT_DTYPE)
        arr = arr.reshape(len(arr) // elementsize, elementsize)
        np.save(output_file, arr)


if __name__ == '__main__':
#    (signal_file, signal_trim_length, catalog_input,
#     output_file, num_parallel) = (
#     '../MXG3.1/adapter-dumps/inventory.h5', 350,
#     '../tables/selected-signal-matches-MXG3.1.txt',
#     'tr.npy', 8)

    (signal_file, signal_trim_length, catalog_input,
     output_file, num_parallel) = sys.argv[1:]

    signal_trim_length = int(signal_trim_length)
    num_parallel = int(num_parallel)

    main(signal_file, signal_trim_length, catalog_input, output_file,
         num_parallel)

