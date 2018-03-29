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
import tarfile
from setup_build import (
    prepare_external_source_dir, build_hdf5, prepare_htslib,
    nanopolish_source_files, htslib_source_files
)

HDF5_URL = 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar.gz'
EIGEN_URL = 'http://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2'
NANOPOLISH_DIR = 'contrib/nanopolish'
HTSLIB_DIR = 'contrib/htslib'
BUILD_DIR = 'build'

prepare_external_source_dir('contrib/eigen', EIGEN_URL, 'eigen', 'bz2')
prepare_external_source_dir('contrib/hdf5', HDF5_URL, 'hdf5', 'gz')
build_hdf5('contrib/hdf5', BUILD_DIR)
prepare_htslib('contrib/htslib')

extmod = Extension('octopus.npinterface',
                   define_macros=[('EIGEN_MPL2_ONLY', 1)],
                   include_dirs=[BUILD_DIR + '/include'] + [
                    os.path.join(NANOPOLISH_DIR, d)
                    for d in """
                        ../htslib ../fast5/include ../eigen src src/hmm
                        src/thirdparty src/thirdparty/scrappie src/common
                        src/alignment src/pore_model eigen""".split()],
                   extra_compile_args=['-fopenmp', '-fsigned-char', '-std=c++11'],
                   library_dirs=[BUILD_DIR + '/lib'],
                   libraries=['z', 'rt', 'dl'], #, 'z', 'rt', 'dl'],
                   extra_link_args=[
                        '-Wl,--strip-all', BUILD_DIR + '/lib/libhdf5.a',
                        '-fopenmp'],
                   sources=["src/npinterface.cc"] + [
                    os.path.join(NANOPOLISH_DIR, d)
                    for d in nanopolish_source_files] + [
                    os.path.join(HTSLIB_DIR, d)
                    for d in htslib_source_files]
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
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=[
        'pomegranate >= 0.9.0',
        'snakemake >= 4.4.0',
        'ont-fast5-api >= 0.4.1',
        'PyYAML >= 3.0',
        'numpy >= 1.14',
    ],
    entry_points={
        'console_scripts': [
            'octopus = octopus.demux_script:main'
        ],
    },
    ext_modules=[extmod],
)
