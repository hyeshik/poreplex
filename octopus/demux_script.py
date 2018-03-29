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

import argparse
import sys
import os
import yaml
from concurrent import futures
from . import __version__
from .signal_analyzer import SignalAnalyzer


def errx(msg):
    print(msg, file=sys.stderr)
    sys.exit(254)


def process_file(inputfile):
    analyzer = SignalAnalyzer(inputfile)

    from octopus.npinterface import get_calibration

    print(inputfile, *get_calibration(inputfile))


def process_batch(inputfiles):
    for f5file in inputfiles:
        process_file(f5file)


def show_banner():
    print("""
\x1b[1mOctopus\x1b[0m version {version} by Hyeshik Chang
- A demultiplexer for nanopore direct RNA sequencing
""".format(version=__version__))


def main_loop(args):
    if not args.quiet:
        show_banner()

    if not os.path.isdir(args.input):
        errx('ERROR: Cannot open the input directory {}.'.format(args.input))

    if not os.path.isdir(args.output):
        try:
            os.makedirs(args.output)
        except:
            errx('ERROR: Failed to create the output directory {}.'.format(args.output))

    if not args.config:
        print(__module__)

    raise SystemExit


    no_parallel = args.parallel <= 1

    with futures.ProcessPoolExecutor(args.parallel) as executor:
        jobs = []
        batch_queue = []

        # Scan FAST5 files and send job chunks to the worker queue.
        def flush():
            if batch_queue:
                if no_parallel:
                    job = process_batch(batch_queue[:])
                else:
                    job = executor.submit(process_batch, batch_queue[:])
                jobs.append(job)
                del batch_queue[:]

        for path, dirs, files in os.walk(args.input):
            for fn in files:
                if fn.lower().endswith('.fast5'):
                    f5path = os.path.join(path, fn)
                    batch_queue.append(f5path)
                    if len(batch_queue) >= args.batch_chunk:
                        flush()
        else:
            flush()

        for job in jobs:
            if job is not None:
                job.result() # TODO: interpret the result.

    print('Done.')

def main():
    parser = argparse.ArgumentParser(
        prog='octopus',
        description='Barcode demultiplexer for nanopore direct RNA sequencing')

    parser.add_argument('-i', '--input', required=True,
                        help='Path to the directory with the input FAST5 files.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory path')
    parser.add_argument('-p', '--parallel', default=1, type=int,
                        help='Number of worker processes (default: 1)')
    parser.add_argument('--batch-chunk', default=32, type=int,
                        help='Number of files in a single batch (default: 32)')
    parser.add_argument('-c', '--config', default='',
                        help='Path to signal processing configuration.')
    parser.add_argument('-q', '--quiet', default=False, action='store_true',
                        help='Suppress non-error messages.')

    args = parser.parse_args(sys.argv[1:])
    main_loop(args)

