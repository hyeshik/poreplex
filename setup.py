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

extra_dependencies_barcoding = ['tensorflow >= 1.8.0', 'Keras >= 2.1.6']
extra_dependencies_live = ['inotify >= 0.2.9']

setup(
    name='octopus',
    packages=['octopus'],
    version='0.1',
    description='Squiggle-level data preprocessor for nanopore direct RNA sequencing',
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
        'h5py >= 2.8.0rc1',
        'ont-fast5-api >= 0.4.1',
        'PyYAML >= 3.0',
        'numpy >= 1.14',
        'pandas >= 0.22.0',
        'scipy >= 1.0',
        'pysam >= 0.14.0',
        'progressbar2 >= 3.37.0',
    ],
    extras_require={
        'barcoding': extra_dependencies_barcoding,
        'live': extra_dependencies_live,
        'full': extra_dependencies_live + extra_dependencies_barcoding,
    }
    entry_points={
        'console_scripts': [
            'octopus = octopus.commandline:__main__'
        ],
    },
)
