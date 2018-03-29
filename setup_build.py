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

from urllib import request
import os
import sys
import tarfile
import re
import shutil
import subprocess as sp

def prepare_external_source_dir(destdir, url, name, compression):
    parentdir = os.path.dirname(destdir)

    tarpath = os.path.join(parentdir, name + '.tar.' + compression)
    if not os.path.exists(tarpath):
        print('Downloading', url)
        request.urlretrieve(url, tarpath)

    if not os.path.exists(os.path.join(destdir, 'CMakeLists.txt')):
        with tarfile.open(tarpath, 'r:' + compression) as archive:
            firstfile = next(iter(archive))
            archive_dirname = firstfile.name.split('/')[0]
            archive.extractall(parentdir)
            os.rename(os.path.join(parentdir, archive_dirname),
                      destdir)

def build_hdf5(sourcedir, destdir):
    if os.path.exists(os.path.join(destdir, 'lib/libhdf5.a')):
        return

    sp.check_call(['env', 'CFLAGS=-fPIC -DPIC -O2',
                   './configure', '--prefix=' + os.path.abspath(destdir),
                   '--enable-threadsafe', '--disable-shared',
                   '--disable-hl'], cwd=sourcedir)
    sp.check_call(['make', '-j4'], cwd=sourcedir)
    sp.check_call(['make', 'install'], cwd=sourcedir)

def prepare_htslib(sourcedir):
    sp.check_call(['make', 'version.h'], cwd=sourcedir)

    # Disable LZMA and BZ support. We don't use them.
    open(os.path.join(sourcedir, 'config.h'), 'w').write("""\
#define HAVE_FSEEKO 1
#define HAVE_DRAND48 1
""")

#    if os.path.exists(os.path.join(destdir, 'lib/libhts.a')):
#        return
#
#    destdir_abs = os.path.abspath(destdir)
#    patched_makefile = os.path.join(destdir_abs, 'Makefile.libhts')
#
#    makefile = open(os.path.join(sourcedir, 'Makefile')).read()
#    open(patched_makefile, 'w').write(
#        re.sub('(\nCFLAGS.*)', lambda m: m.group() + ' -fPIC -DPIC',
#               makefile))
#    sp.check_call(['make', '-f', patched_makefile, '-j4'], cwd=sourcedir)


nanopolish_source_files = """
    src/nanopolish_variant_db.cpp
    src/nanopolish_getmodel.cpp
    src/nanopolish_haplotype.cpp
    src/nanopolish_call_methylation.cpp
    src/nanopolish_squiggle_read.cpp
    src/nanopolish_methyltrain.cpp
    src/nanopolish_raw_loader.cpp
    src/nanopolish_index.cpp
    src/nanopolish_scorereads.cpp
    src/training_core.cpp
    src/nanopolish_phase_reads.cpp
    src/nanopolish_extract.cpp
    src/nanopolish_read_db.cpp
    src/nanopolish_call_variants.cpp
    src/nanopolish_train_poremodel_from_basecalls.cpp
    src/hmm/nanopolish_transition_parameters.cpp
    src/hmm/nanopolish_duration_model.cpp
    src/hmm/nanopolish_profile_hmm_r9.cpp
    src/hmm/nanopolish_profile_hmm.cpp
    src/hmm/nanopolish_profile_hmm_r7.cpp
    src/common/fs_support.cpp
    src/common/nanopolish_bam_processor.cpp
    src/common/nanopolish_klcs.cpp
    src/common/nanopolish_alphabet.cpp
    src/common/logsum.cpp
    src/common/nanopolish_iupac.cpp
    src/common/nanopolish_common.cpp
    src/common/nanopolish_variant.cpp
    src/common/nanopolish_bam_utils.cpp
    src/common/nanopolish_fast5_io.cpp
    src/alignment/nanopolish_eventalign.cpp
    src/alignment/nanopolish_anchor.cpp
    src/alignment/nanopolish_alignment_db.cpp
    src/pore_model/nanopolish_model_names.cpp
    src/pore_model/nanopolish_pore_model_set.cpp
    src/pore_model/nanopolish_poremodel.cpp
    src/thirdparty/stdaln.c
    src/thirdparty/scrappie/event_detection.c
    src/thirdparty/scrappie/scrappie_common.c
    src/thirdparty/scrappie/util.c
""".split()

htslib_source_files = """
    hts.c hfile.c hfile_net.c bgzf.c
    thread_pool.c cram/pooled_alloc.c faidx.c kstring.c
    sam.c cram/cram_io.c
    cram/string_alloc.c cram/sam_header.c
    cram/open_trace_file.c
    cram/cram_index.c
    cram/cram_encode.c
    multipart.c textutils.c md5.c cram/mFILE.c
    cram/cram_stats.c cram/cram_decode.c
    knetfile.c
    cram/cram_samtools.c
    cram/files.c
    hts_os.c
    cram/cram_codecs.c
    cram/rANS_static.c
""".split()
