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

from pysam import BGZFile, faidx
from glob import glob
from functools import partial
from collections import defaultdict
from shutil import copyfileobj
import h5py
import numpy as np
import logging
import os
from errno import EXDEV, EEXIST


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
            if entry.get('fastq') is not None:
                seq, qual, adapter_length = entry['fastq']
                if adapter_length > 0:
                    seq = seq[:-adapter_length]
                    qual = qual[:-adapter_length]

                formatted = '@{}\n{}\n+\n{}\n'.format(entry['read_id'], seq, qual)
                self.streams[entry['label']].write(formatted.encode('ascii'))


def link_fast5_files(config, results):
    indir = config['inputdir']
    outdir = config['outputdir']
    symlinkfirst = config['fast5_always_symlink']
    labelmap = config['output_names']
    blacklist_hardlinks = set()

    for entry in results:
        if 'label' not in entry: # Error before opening the FAST5
            continue

        original_fast5 = os.path.join(indir, entry['filename'])
        link_path = os.path.join(outdir, 'fast5', labelmap[entry['label']],
                                 entry['filename'])

        original_dir = os.path.dirname(original_fast5)
        link_dir = os.path.dirname(link_path)

        if not os.path.isdir(link_dir):
            try:
                os.makedirs(link_dir)
            except FileExistsError:
                pass

        for retry in [0, 1]:
            if not symlinkfirst and (original_dir, link_dir) not in blacklist_hardlinks:
                try:
                    os.link(original_fast5, link_path)
                except OSError as exc:
                    if exc.errno == EEXIST:
                        os.unlink(link_path)
                        continue
                    elif exc.errno != EXDEV:
                        raise
                    blacklist_hardlinks.add((original_dir, link_dir))
                else:
                    break

            try:
                os.symlink(os.path.abspath(original_fast5), link_path)
            except OSError as exc:
                if exc.errno == EEXIST:
                    os.unlink(link_path)
                    continue
                else:
                    raise
            else:
                break


class SequencingSummaryWriter:

    SUMMARY_OUTPUT_FIELDS = [
        'filename', 'read_id', 'run_id', 'channel', 'start_time',
        'duration', 'num_events', 'sequence_length', 'mean_qscore',
        'sample_id', 'status', 'label',
    ]

    def __init__(self, outputdir, labelmapping):
        self.file = open(os.path.join(outputdir, 'sequencing_summary.txt'), 'w')
        self.labelmapping = labelmapping
        print(*self.SUMMARY_OUTPUT_FIELDS, sep='\t', file=self.file)

    def close(self):
        self.file.close()

    def write_results(self, results):
        for entry in results:
            if 'label' in entry:
                print(*[self.labelmapping.get(entry[f], entry[f])
                        if f == 'label' else entry[f]
                        for f in self.SUMMARY_OUTPUT_FIELDS],
                      file=self.file, sep='\t')


class NanopolishReadDBWriter:

    def __init__(self, outputdir, labelmapping):
        self.labelmapping = labelmapping
        self.outputdir = os.path.join(outputdir, 'nanopolish')
        self.seqfiles, self.dbfiles = self.open_streams()

    def open_streams(self):
        seqfiles, dbfiles = {}, {}
        for groupid, name in self.labelmapping.items():
            filepath = os.path.join(self.outputdir, name + '.fasta')
            seqfiles[groupid] = open(filepath, 'w')
            dbfiles[groupid] = open(filepath + '.index.readdb', 'w')
        return seqfiles, dbfiles

    def close(self):
        for groupid, f in list(self.seqfiles.items()):
            f.close()
            del self.seqfiles[groupid]

        for groupid, f in list(self.dbfiles.items()):
            f.close()
            del self.dbfiles[groupid]

        # Create bgzipped-fasta and indices
        for groupid, name in self.labelmapping.items():
            inputfile = os.path.join(self.outputdir, name + '.fasta')
            if os.path.getsize(inputfile) > 0:
                bgzippedfile = inputfile + '.index'
                copyfileobj(open(inputfile, 'rb'), BGZFile(bgzippedfile, 'w'))
                faidx(bgzippedfile)

    def write_sequences(self, procresult):
        for entry in procresult:
            if entry.get('fastq') is not None:
                formatted = '>{}\n{}\n'.format(entry['read_id'], entry['fastq'][0])
                self.seqfiles[entry['label']].write(formatted)

                fast5_relpath = os.path.join('..', 'fast5', self.labelmapping[entry['label']],
                                             entry['filename'])
                formatted = '{}\t{}\n'.format(entry['read_id'], fast5_relpath)
                self.dbfiles[entry['label']].write(formatted)


class FinalSummaryTracker:

    REPORTING_ORDER = ['pass', 'artifact', 'fail', 'file_error']
    FRIENDLY_LABELS = {
        0: 'Barcoded sample 1 (BC1)',
        1: 'Barcoded sample 2 (BC2)',
        2: 'Barcoded sample 3 (BC3)',
        3: 'Barcoded sample 4 (BC4)',
        'pass': 'Successfully processed',
        'fail': 'Processing failed',
        'artifact': 'Possibly artifact',
        'file_error': 'Failed to open',
    }
    PASS_LABEL_WITH_BARCODES = 'Barcode undetermined'
    FRIENDLY_STATUS = {
        'fail': {
            'adapter_not_detected': "3' Adapter could not be located",
            'not_basecalled': 'No albacore basecall data found',
            'too_few_events': 'Signal is too short',
        },
        'artifact': {
            'unsplit_read': 'Two or more molecules found in a single read',
        },
        'file_error': {
            'disappeared': 'File is moved to other location',
            'unknown_error': 'File could not be opened due to unknown error',
        },
    }

    def __init__(self, labelmapping):
        self.labelmapping = labelmapping
        self.counts = defaultdict(int)

        barcodenums = [k for k in self.labelmapping.keys() if not isinstance(k, str)]
        self.reporting_order = (
            sorted(barcodenums) + self.REPORTING_ORDER
            if barcodenums else self.REPORTING_ORDER)

        if barcodenums:
            self.FRIENDLY_LABELS['pass'] = self.PASS_LABEL_WITH_BARCODES

    def feed_results(self, results):
        for entry in results:
            self.counts[entry.get('label', 'file_error'), entry['status']] += 1

    def print_results(self, file):
        if hasattr(file, 'write'):
            _ = partial(print, sep='\t', file=file)
        else:
            logger = logging.getLogger('poreplex')
            _ = lambda *args: logger.error(' '.join(map(str, args)))

        _("== Result Summary ==")
        for label in self.reporting_order:
            matching_subcounts = {s: cnt for (l, s), cnt in self.counts.items() if l == label}
            subtotal = sum(matching_subcounts.values())

            friendlylabel = self.FRIENDLY_LABELS[label]
            _(" * {}:{}".format(friendlylabel, '\t' if len(friendlylabel) < 20 else ''), subtotal)
            if label in self.FRIENDLY_STATUS:
                for status, cnt in matching_subcounts.items():
                    _("    - {}:".format(self.FRIENDLY_STATUS[label][status]), cnt)

        _('')

def create_links_rebalanced(desth5, group, infiles):
    desth5.require_group(group)

    for datafile in infiles:
        basename = os.path.basename(datafile)
        with h5py.File(datafile, 'r') as d5:
            for batchid, subgrp in d5[group].items():
                for readid in subgrp.keys():
                    dumpgroup = get_read_id_dump_group(readid)
                    gobj = desth5.require_group(group + '/' + dumpgroup)
                    try:
                        gobj[readid] = h5py.ExternalLink(basename,
                            '{}/{}/{}'.format(group, batchid, readid))
                    except RuntimeError:
                        if readid not in gobj:
                            raise

def create_adapter_dumps_inventory(destfile, filepattern):
    with h5py.File(destfile, 'w') as ivt:
        # Merge items under catalog/adapter group.
        ivt.require_group('catalog')

        fragments = []
        for datafile in glob(filepattern):
            with h5py.File(datafile, 'r') as d5:
                for batchid, tbl in d5['catalog/adapter'].items():
                    fragments.append(tbl[:])

        fulltbl = np.hstack(fragments)
        fulltbl.sort(order='read_id')
        ivt['catalog/adapter'] = fulltbl

        # Create hashed groups of external links for the signal datasets.
        create_links_rebalanced(ivt, 'adapter', glob(filepattern))

def create_events_inventory(destfile, filepattern):
    with h5py.File(destfile, 'w') as ivt:
        create_links_rebalanced(ivt, 'basecalled_events', glob(filepattern))

def get_read_id_dump_group(read_id, grplength=3):
    #hash = unpack('I', sha1(read_id.encode('ascii')).digest()[:4])[0] % hashsize
    #return format(hash, numformat)
    return read_id[:grplength]
