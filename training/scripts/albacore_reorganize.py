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

from ont_fast5_api.fast5_file import Fast5File
from ont_fast5_api.analysis_tools.basecall_1d import Basecall1DTools
from concurrent import futures
from collections import defaultdict
import pandas as pd
import gzip
import os

META_FILES = """
configuration.cfg
pipeline.log
sequencing_telemetry.js""".split()


def get_fast5_readid(filename):
    with Fast5File(filename, 'r') as f5:
        readinfo = f5.status.read_info[0]

        try:
            with Basecall1DTools(f5, 'r') as bc:
                _, sequence, qscore = bc.get_called_sequence('template')
        except KeyError:
            sequence = qscore = ''

        return (f5.get_tracking_id()['run_id'],
                readinfo.read_id,
                f5.get_channel_info()['channel_number'],
                readinfo.read_number,
                readinfo.duration,
                sequence, qscore)

def scan_directory(dirpath, files):
    results = []

    print('Processing', dirpath)

    for fn in files:
        if fn.endswith('.fast5'):
            f5path = os.path.join(dirpath, fn)
            datafields = get_fast5_readid(f5path)
            results.append((f5path,) + datafields)

    return results

def merge_sequencing_summary(outdir, srcdir):
    inpath = os.path.join(srcdir, 'sequencing_summary.txt')
    statspath = os.path.join(outdir, 'sequences.txt.gz')

    seqsummary = pd.read_table(inpath, low_memory=False)
    seqstats = pd.read_table(statspath, compression='gzip', low_memory=False)
    seqmg = pd.merge(seqsummary, seqstats, how='right', left_on='read_id',
                     right_on='read_id', suffixes=['_x', ''])

    outpath = os.path.join(outdir, 'sequencing_summary.txt.gz')
    seqmg[seqsummary.columns.tolist() + ['read_no']
         ].to_csv(outpath, sep='\t', compression='gzip', index=False)

    seqmissing = set(seqsummary['read_id'].tolist()) - set(seqmg['read_id'].tolist())
    seqmissing = seqsummary[seqsummary['read_id'].isin(seqmissing)]
    missingout = os.path.join(outdir, 'sequencing_summary_missing.txt.gz')
    seqmissing.to_csv(missingout, sep='\t', compression='gzip', index=False)

def copy_meta_files(outdir, srcdir):
    for mf in META_FILES:
        inpath = os.path.join(srcdir, mf)
        outpath = os.path.join(outdir, mf + '.gz')
        gzip.open(outpath, 'wt').write(open(inpath).read())

def main(outdir, srcdir, nproc):
    runidmap = defaultdict(lambda: len(runidmap))

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    fastq_out = gzip.open(os.path.join(outdir, 'sequences.fastq.gz'), 'wt')
    fasta_out = gzip.open(os.path.join(outdir, 'sequences.fasta.gz'), 'wt')
    catalog_out = gzip.open(os.path.join(outdir, 'sequences.txt.gz'), 'wt')
    print('filename', 'run_id', 'read_id', 'channel', 'read_no', 'duration',
          'sequence_length', sep='\t', file=catalog_out)

    try:
        with futures.ProcessPoolExecutor(nproc) as executor:
            jobs = []

            for dirpath, dirs, files in os.walk(srcdir):
                job = executor.submit(scan_directory, dirpath, files)
                jobs.append(job)

            for job in jobs:
                for origfile, runid, readid, ch, readno, duration, seq, qscore in job.result():
                    runno = runidmap[runid]
                    newf5dir = os.path.join(outdir, 'fast5/{}/ch{}'.format(runno, ch))
                    newf5path = '{}/r{}.fast5'.format(newf5dir, readno)
                    if not os.path.isdir(newf5dir):
                        os.makedirs(newf5dir)
                    os.link(origfile, newf5path)

                    print(newf5path, runid, readid, ch, readno, duration, len(seq),
                          sep='\t', file=catalog_out)
                    if len(seq) > 0:
                        print('>{}\n{}'.format(readid, seq), file=fasta_out)
                        print('@{}\n{}\n+\n{}'.format(readid, seq, qscore), file=fastq_out)
    finally:
        fastq_out.close()
        fasta_out.close()
        catalog_out.close()

    merge_sequencing_summary(outdir, srcdir)
    copy_meta_files(outdir, srcdir)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Tidy up the albacore output directory.')
    parser.add_argument('-i', '--input', required=True,
			help='path to albacore output files')
    parser.add_argument('-o', '--output', required=True,
			help='path to save reorganized files')
    parser.add_argument('-p', '--parallel', type=int, default=1,
			help='number of parallel processes (default: 1)')

    args = parser.parse_args()
    return (args.input, args.output, args.parallel)

if __name__ == '__main__':
    indir, outdir, nproc = parse_arguments()
    main(outdir, indir, nproc)

