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

from pysam import BGZFile
from glob import glob
import h5py
import os


class FASTQWriter:

    def __init__(self, outputdir, names):
        self.outputdir = outputdir
        self.names = names

        self.open_streams()

    def open_streams(self):
        self.streams = {
            int_name: BGZFile(self.get_output_path(name), 'w')
            for int_name, name in self.names.items()}

    def close(self):
        for stream in self.streams.values():
            stream.close()

    def get_output_path(self, name):
        return os.path.join(self.outputdir, 'fastq', name + '.fastq.gz')

    def write_sequences(self, procresult):
        for entry in procresult:
            if entry['fastq'] is not None:
                formatted = ''.join('@{}\n{}\n+\n{}\n'.format(entry['read_id'], *entry['fastq']))
                self.streams[entry['label']].write(formatted.encode('ascii'))


class SequencingSummaryWriter:

    SUMMARY_OUTPUT_FIELDS = [
        'filename', 'read_id', 'run_id', 'channel', 'start_time',
        'duration', 'num_events', 'sequence_length', 'mean_qscore',
        'sample_id', 'error', 'label',
    ]

    def __init__(self, outputdir):
        self.file = open(os.path.join(outputdir, 'sequencing_summary.txt'), 'w')
        print(*self.SUMMARY_OUTPUT_FIELDS, sep='\t', file=self.file)

    def close(self):
        self.file.close()

    def write_results(self, results):
        for entry in results:
            print(*[entry[f] for f in self.SUMMARY_OUTPUT_FIELDS],
                  file=self.file, sep='\t')

def create_inventory_hdf5(destfile, filepattern, groupnames):
    def link_group_items(groupname):
        destgrp = ivt.require_group(groupname)

        for datafile in glob(filepattern):
            basename = os.path.basename(datafile)
            with h5py.File(datafile, 'r') as d5:
                for k in d5[groupname].keys():
                    destgrp[k] = h5py.ExternalLink(basename, groupname + '/' + k)

    with h5py.File(destfile, 'w') as ivt:
        for groupname in groupnames:
            link_group_items(groupname)
