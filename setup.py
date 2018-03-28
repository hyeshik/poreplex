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

from setuptools import setup
from distutils.core import Extension
from urllib import request
import os
import sys
import pkgconfig
import tarfile

EIGEN_URL = 'http://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2'
NANOPOLISH_DIR = 'contrib/nanopolish'

# If nanopolish module is linked with the bundled libraries, the shared
# libraries brought by h5py and pysam can conflict against those.
HTS_CFLAGS = pkgconfig.cflags('htslib')
HTS_LIBS = pkgconfig.libs('htslib')
HDF5_CFLAGS = pkgconfig.cflags('hdf5')
HDF5_LIBS = pkgconfig.libs('hdf5')

if not HTS_LIBS:
    print('ERROR: pkg-config could not detect htslib.', file=sys.stderr)
if not HDF5_LIBS:
    print('ERROR: pkg-config could not detect libhdf5.', file=sys.stderr)
if not HTS_LIBS or not HDF5_LIBS:
    sys.exit(254)

def split_dir_switches(fullstr, *prefixes):
    rest = []
    prefixed = [[] for _ in prefixes]

    for sw in fullstr.split():
        for pstack, prefix in zip(prefixed, prefixes):
            if sw.startswith(prefix):
                pstack.append(sw[len(prefix):])
                break
        else:
            rest.append(sw)

    return prefixed + [rest]

HTS_INCLUDE_DIRS, HTS_EXTRA_CFLAGS = split_dir_switches(HTS_CFLAGS, '-I')
HTS_LIBRARY_DIRS, HTS_SHARED_LIBS, HTS_EXTRA_LIBS = (
    split_dir_switches(HTS_LIBS, '-L', '-l'))
HDF5_INCLUDE_DIRS, HDF5_EXTRA_CFLAGS = split_dir_switches(HDF5_CFLAGS, '-I')
HDF5_LIBRARY_DIRS, HDF5_SHARED_LIBS, HDF5_EXTRA_LIBS = (
    split_dir_switches(HDF5_LIBS, '-L', '-l'))

def prepare_eigen_source_dir(destdir):
    parentdir = os.path.dirname(destdir)

    tarpath = os.path.join(parentdir, 'eigen.tar.bz2')
    if not os.path.exists(tarpath):
        print('Downloading', EIGEN_URL)
        request.urlretrieve(EIGEN_URL, tarpath)

    if not os.path.exists(os.path.join(destdir, 'INSTALL')):
        with tarfile.open(tarpath, 'r:bz2') as archive:
            firstfile = next(iter(archive))
            archive_dirname = firstfile.name.split('/')[0]
            archive.extractall(parentdir)
            os.rename(os.path.join(parentdir, archive_dirname),
                      destdir)

prepare_eigen_source_dir('contrib/eigen')

extmod = Extension('octopus.npinterface',
                   define_macros=[('EIGEN_MPL2_ONLY', 1)],
                   include_dirs=HTS_INCLUDE_DIRS + HDF5_INCLUDE_DIRS + [
                    os.path.join(NANOPOLISH_DIR, d)
                    for d in """
                        ../fast5/include ../eigen src src/hmm
                        src/thirdparty src/thirdparty/scrappie src/common
                        src/alignment src/pore_model eigen""".split()],
                   extra_compile_args=HTS_EXTRA_CFLAGS + HDF5_EXTRA_CFLAGS + [
                        '-fopenmp', '-fsigned-char', '-std=c++11'],
                   library_dirs=HTS_LIBRARY_DIRS + HDF5_LIBRARY_DIRS,
                   libraries=HTS_SHARED_LIBS + HDF5_SHARED_LIBS, #, 'z', 'rt', 'dl'],
                   extra_link_args=HTS_EXTRA_LIBS + HDF5_EXTRA_LIBS + [
                        '-fopenmp'],
                   sources=["src/npinterface.cc"] + [
                    os.path.join(NANOPOLISH_DIR, d)
                    for d in """
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
                        """.split()]
                  )

setup(
    name='octopus',
    packages=['octopus'],
    version='0.1',
    description='Barcode demultiplexer for nanopore direct RNA sequencing',
    author='Hyeshik Chang',
    author_email='hyeshik@snu.ac.kr',
    url='https://github.com/hyeshik/octopus',
    download_url='https://github.com/hyeshik/octopus/archive/octopus-0.1.tar.gz',
    keywords=[
        'nanopore',
        'direct RNA sequencing',
        'barcode',
        'demultiplexing'
    ],
    license='MIT',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=[
        'pomegranate >= 0.9.0',
        'snakemake >= 4.4.0',
        'ont-fast5-api >= 0.4.1',
    ],
    entry_points={
        'console_scripts': [
            'octopus = octopus.demux_script:main'
        ],
    },
    ext_modules=[extmod],
)
