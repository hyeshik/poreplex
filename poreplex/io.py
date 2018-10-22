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
from threading import Lock
import h5py
import numpy as np
import pandas as pd
import logging
import os
from errno import EXDEV, EEXIST
from .utils import ensure_dir_exists


class FASTQWriter:

    def __init__(self, output_dir, output_layout):
        self.output_dir = output_dir
        self.output_layout = output_layout
        self.lock = Lock()

        self.open_streams()

    def open_streams(self):
        self.streams = {
            int_name: BGZFile(self.get_output_path(name), 'w')
            for int_name, name in self.output_layout.items()}

    def close(self):
        for stream in self.streams.values():
            stream.close()

    def get_output_path(self, name):
        output_path = os.path.join(self.output_dir, 'fastq', name + '.fastq.gz')
        ensure_dir_exists(output_path)
        return output_path

    def write_sequences(self, procresult):
        with self.lock:
            for entry in procresult:
                if entry.get('sequence') is not None:
                    seq, qual, adapter_length = entry['sequence']
                    if adapter_length > 0:
                        seq = seq[:-adapter_length]
                        qual = qual[:-adapter_length]

                    output_name = entry['label'], entry.get('barcode')
                    formatted = '@{}\n{}\n+\n{}\n'.format(entry['read_id'], seq, qual)
                    self.streams[output_name].write(formatted.encode('ascii'))


def link_fast5_files(config, results):
    indir = config['inputdir']
    outdir = config['outputdir']
    symlinkfirst = config['fast5_always_symlink']
    output_layout = config['output_layout']
    blacklist_hardlinks = set()

    for entry in results:
        if 'label' not in entry: # Error before opening the FAST5
            continue

        original_fast5 = os.path.join(indir, entry['filename'])
        link_path = os.path.join(outdir, 'fast5',
                                 output_layout[entry['label'], entry.get('barcode')],
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

    def __init__(self, config, output_dir, label_mapping, barcode_mapping):
        self.file = open(os.path.join(output_dir, 'sequencing_summary.txt'), 'w')
        self.lock = Lock()
        self.label_mapping = label_mapping

        self.output_fields = self.SUMMARY_OUTPUT_FIELDS[:]

        if config['barcoding']:
            self.barcode_mapping = barcode_mapping
            self.output_fields.append('barcode')
        else:
            self.barcode_mapping = None

        if config['measure_polya']:
            self.polya_enabled = True
            self.output_fields.append('polya_dwell')
        else:
            self.polya_enabled = False

        if config['fast5_output']:
            if config['barcoding']:
                self.format_filename = (lambda entry:
                    os.path.join('fast5', entry['label'],
                                 self.barcode_mapping[entry.get('barcode')],
                                 entry['filename']))
            else:
                self.format_filename = (lambda entry:
                    os.path.join('fast5', entry['label'], entry['filename']))
        else:
            self.format_filename = lambda entry: entry['filename']

        print(*self.output_fields, sep='\t', file=self.file)

    def close(self):
        self.file.close()

    def write_results(self, results):
        with self.lock:
            for entry in results:
                if 'label' in entry:
                    output_entry = entry.copy()
                    output_entry['label'] = self.label_mapping[entry['label']]
                    output_entry['filename'] = self.format_filename(output_entry)
                    if self.barcode_mapping is not None:
                        output_entry['barcode'] = self.barcode_mapping[entry.get('barcode')]
                    if self.polya_enabled:
                        output_entry['polya_dwell'] = (
                            format(entry['polya']['dwell_time'], '.4f')
                            if 'polya' in entry else '')

                    print(*[output_entry[f] for f in self.output_fields],
                          file=self.file, sep='\t')


class NanopolishReadDBWriter:

    def __init__(self, output_dir, output_layout):
        self.output_layout = output_layout
        self.output_dir = os.path.join(output_dir, 'nanopolish')
        self.seqfiles, self.dbfiles = self.open_streams()
        self.lock = Lock()

    def open_streams(self):
        seqfiles, dbfiles = {}, {}
        for groupid, name in self.output_layout.items():
            filepath = os.path.join(self.output_dir, name + '.fasta')
            ensure_dir_exists(filepath)
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
        for groupid, name in self.output_layout.items():
            inputfile = os.path.join(self.output_dir, name + '.fasta')
            if os.path.getsize(inputfile) > 0:
                bgzippedfile = inputfile + '.index'
                copyfileobj(open(inputfile, 'rb'), BGZFile(bgzippedfile, 'w'))
                faidx(bgzippedfile)

    def write_sequences(self, procresult):
        with self.lock:
            for entry in procresult:
                if entry.get('sequence') is not None:
                    mappingkey = entry['label'], entry.get('barcode')

                    formatted = '>{}\n{}\n'.format(entry['read_id'], entry['sequence'][0])
                    self.seqfiles[mappingkey].write(formatted)

                    fast5_relpath = os.path.join('fast5', self.output_layout[mappingkey],
                                                 entry['filename'])
                    formatted = '{}\t{}\n'.format(entry['read_id'], fast5_relpath)
                    self.dbfiles[mappingkey].write(formatted)


class FinalSummaryTracker:

    REPORTING_ORDER = ['pass', 'artifact', 'fail']
    FRIENDLY_LABELS = {
        'pass': 'Successfully processed',
        'fail': 'Processing failed',
        'artifact': 'Possibly artifact',
    }
    BARCODE_FRIENDLY_NAME = 'Barcoded sample {num} (BC{num})'
    FRIENDLY_STATUS = {
        'fail': {
            'scaler_signal_too_short': 'Signal is too short',
            'sequence_too_short': 'Sequence is too short',
            'irregular_fast5': 'Invalid FAST5 format',
            'basecall_table_incomplete': 'Basecall table does not match',
            'adapter_not_detected': "3' Adapter could not be located",
            'not_basecalled': 'No albacore basecall data found',
            'scaling_qc_fail': 'Signal scaling QC failed',
            'disappeared': 'File is moved to other location',
            'unknown_error': 'File could not be opened due to unknown error',
        },
        'artifact': {
            'unsplit_read': 'Two or more molecules found within a read',
        },
    }

    LABEL_FORMAT = '{:49s} '
    LABEL_BULLET = ' - '
    MINIMUM_COLUMN_WIDTH = 3

    def __init__(self, label_names, barcode_names):
        self.label_names = label_names
        self.barcode_names = barcode_names
        self.counts = defaultdict(int)
        self.label_reporting_order = self.REPORTING_ORDER
        self.barcode_reporting_order = sorted([
            n for n in barcode_names.keys() if n is not None]) + [None]

    def feed_results(self, results):
        for entry in results:
            self.counts[entry.get('label', 'fail'),
                        entry.get('barcode', None),
                        entry['status']] += 1

    def print_results(self, file):
        if hasattr(file, 'write'):
            _ = partial(print, sep='\t', file=file)
        else:
            logger = logging.getLogger('poreplex')
            _ = lambda *args: logger.error(' '.join(map(str, args)))

        _("==== Result Summary ====")
        longest_count_length = len(format(max(self.counts.values()), 'd'))
        column_width = max(self.MINIMUM_COLUMN_WIDTH, longest_count_length)
        column_title_format = '{{:{}s}} '.format(column_width)
        column_numeric_format = '{{:{}d}} '.format(column_width)

        # Show the header
        if len(self.barcode_names) > 1:
            label_fields = [self.LABEL_FORMAT.format('')] + [
                column_title_format.format(self.barcode_names[bc])
                for bc in self.barcode_reporting_order
            ]
            _(''.join(label_fields))

        # Show the table content
        tbl = pd.DataFrame([(k[0], -1 if k[1] is None else k[1], k[2], v)
                            for k, v in self.counts.items()],
                           columns=['label', 'barcode', 'status', 'count'])
        tbl['label_key'] = tbl['label'].apply(self.label_reporting_order.index)
        ordered = (tbl.sort_values(by=['label_key', 'count'], ascending=[True, False])
                      .groupby(by=['label', 'status'], sort=False))
        current_label = None
        for lk, grp in ordered:
            linelabel = None

            # Handle label changes
            if current_label is None or current_label != lk[0]:
                current_label = lk[0]
                if current_label in self.FRIENDLY_STATUS:
                    _(self.LABEL_FORMAT.format(self.FRIENDLY_LABELS[current_label]))
                else:
                    linelabel = self.FRIENDLY_LABELS[current_label]

            # Set the line title for the rows with detailed status explanations
            if linelabel is None:
                linelabel = self.LABEL_BULLET + self.FRIENDLY_STATUS[current_label][lk[1]]

            readcount_by_barcode = grp.set_index('barcode')['count'].to_dict()
            readcounts = [readcount_by_barcode.get(bc if bc is not None else -1, 0)
                          for bc in self.barcode_reporting_order]

            # Print the read counts
            _(self.LABEL_FORMAT.format(linelabel) +
              ''.join(column_numeric_format.format(cnt) for cnt in readcounts))

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
