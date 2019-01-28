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

import mappy
from pysam import FUNMAP, FREVERSE, FSECONDARY, FSUPPLEMENTARY
from pysam import AlignmentFile, AlignedSegment
from struct import pack, unpack, calcsize
from collections import defaultdict
from threading import Lock
from .utils import ensure_dir_exists

MM_IDX_MAGIC = b"MMI\2"


def check_minimap2_index(filename):
    with open(filename, 'rb') as idxf:
        magic = idxf.read(4)
        if magic != MM_IDX_MAGIC:
            raise Exception('File magic is not found from ' + filename)


class BAMWriter:

    def __init__(self, output, indexed_sequence_list, index_options):
        header = self.build_header(indexed_sequence_list, index_options)
        ensure_dir_exists(output)
        self.writer = AlignmentFile(output, 'wb', header=header)
        self.lock = Lock()

    def __del__(self):
        self.close()

    def close(self):
        if hasattr(self, 'writer'):
            self.writer.close()
            del self.writer

    def build_header(self, indexed_sequence_list, index_options):
        return {'SQ': indexed_sequence_list,
                'PG': [{'ID': 'minimap2', 'PN': 'minimap2',
                        'CL': index_options, 'DS': 'minimap2 invoked by poreplex'}]}

    def write(self, fields):
        line = '\t'.join(map(str, fields))
        segment = AlignedSegment.fromstring(line, self.writer.header)
        with self.lock:
            self.writer.write(segment)


class AlignmentWriter:

    def __init__(self, indexfile, output, output_layout):
        self.aligner = mappy.Aligner(indexfile)
        if not self.aligner:
            raise Exception('Could not open minimap2 index {}.'.format(indexfile))
        self.writers = self.open_writers(indexfile, output, output_layout)

    def open_writers(self, indexfile, output, output_layout):
        indexed_sequences, index_options = list(self.get_indexed_sequence_list(indexfile))
        return {muxid: BAMWriter(output.format(name), indexed_sequences, index_options)
                for muxid, name in output_layout.items()}

    def close(self):
        for muxid, writer in self.writers.items():
            writer.close()
        self.writers.clear()

    def __del__(self):
        self.close()

    def get_indexed_sequence_list(self, indexfile):
        seqlist = []

        with open(indexfile, 'rb') as idxf:
            magic = idxf.read(4)
            if magic != MM_IDX_MAGIC:
                raise Exception('File magic is not found from ' + filename)

            header_format = '<IIIII'
            header_size = calcsize(header_format)
            header = idxf.read(header_size)
            if len(header) != header_size:
                raise Exception('Unexpected end of file during reading a header: ' + filename)

            w, k, b, n_seq, flag = unpack(header_format, header)
            index_options = 'minimap2 -w {} -k {}'.format(w, k)
            for i in range(n_seq):
                namlen = idxf.read(1)[0]
                name_seqlen = idxf.read(namlen + 4)
                name = name_seqlen[:-4].decode()
                seqlen = unpack('<I', name_seqlen[-4:])[0]
                seqlist.append({'LN': seqlen, 'SN': name})

        return seqlist, index_options

    def map(self, name, seq, qual):
        seq = seq.replace('U', 'T')
        seqmaps = list(self.aligner.map(seq))
        if not seqmaps:
            yield (name, int(FUNMAP), '*', 0, 0, '*', '*', 0, 0, seq, qual)
            return

        for i, h in enumerate(seqmaps):
            if i > 0:
                flag = int(FSECONDARY)
            elif not h.is_primary:
                flag = int(FSUPPLEMENTARY)
            else:
                flag = 0

            leftclip = '{}S'.format(h.q_st) if h.q_st > 0 else ''
            rightclip = '{}S'.format(len(seq) - h.q_en) if h.q_en < len(seq) else ''

            if h.strand > 0:
                seq_f = seq
                qual_f = qual
            else:
                seq_f = mappy.revcomp(seq)
                qual_f = qual[::-1]
                leftclip, rightclip = rightclip, leftclip
                flag |= FREVERSE

            fullcigar = leftclip + h.cigar_str + rightclip

            yield (name, flag, h.ctg, h.r_st + 1, h.mapq, fullcigar, '*',
                   0, 0, seq_f, qual_f, 'NM:i:{}'.format(h.NM))

    def map_and_write(self, streamid, name, seq, qual, adapter_length):
        writer = self.writers[streamid]
        mapped_seqname = None
        if adapter_length > 0:
            seq = seq[:-adapter_length]
            qual = qual[:-adapter_length]

        for row in self.map(name, seq, qual):
            if mapped_seqname is None:
                mapped_seqname = row[2]
            writer.write(row)
        return mapped_seqname

    def process(self, results):
        mapped_seqs = defaultdict(list)
        failed_counts = defaultdict(int)
        unmapped_counts = defaultdict(int)

        for result in results:
            barcode = result.get('barcode')
            streamid = result.get('label', 'fail'), barcode

            if result.get('sequence') is None or 'read_id' not in result:
                failed_counts[barcode] += 1
            else:
                mapped = self.map_and_write(streamid, result['read_id'], *result['sequence'])

                if mapped == '*':
                    unmapped_counts[barcode] += 1
                else:
                    mapped_seqs[barcode].append(mapped)

        return {'mapped': mapped_seqs, 'failed': failed_counts, 'unmapped': unmapped_counts}

