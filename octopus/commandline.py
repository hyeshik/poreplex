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
import time
import yaml
import shutil
from functools import partial
from . import *
from .pipeline import ProcessingSession
from .utils import *


def show_banner():
    print("""
\x1b[1mOctopus\x1b[0m version {version} by Hyeshik Chang
- A demultiplexer for nanopore direct RNA sequencing
""".format(version=__version__))


def load_config(args):
    presets_dir = os.path.join(os.path.dirname(__file__), 'presets')
    if not args.config:
        config_path = os.path.join(presets_dir, 'rna-r941.cfg')
    elif os.path.isfile(args.config):
        config_path = args.config
    elif os.path.isfile(os.path.join(presets_dir, args.config + '.cfg')):
        config_path = os.path.join(presets_dir, args.config + '.cfg')
    else:
        errx('ERROR: Cannot find a configuration in {}.'.format(args.config))

    config = yaml.load(open(config_path))
    kmer_models_dir = os.path.join(os.path.dirname(__file__), 'kmer_models')
    if not os.path.isabs(config['kmer_model']):
        config['kmer_model'] = os.path.join(kmer_models_dir, config['kmer_model'])

    return config


def create_output_directories(outputdir, config):
    existing = os.listdir(outputdir)
    if existing:
        while config['interactive']:
            answer = input('Output directory {} is not empty. Clear it? (y/N) '
                            .format(outputdir))
            answer = answer.lower()[:1]
            if answer in ('', 'n'):
                sys.exit(1)
            elif answer == 'y':
                print()
                break

        for ent in existing:
            fpath = os.path.join(outputdir, ent)
            if os.path.isdir(fpath):
                shutil.rmtree(fpath)
            else:
                os.unlink(fpath)

    subdirs = ['fastq', 'tmp']

    if config['fast5_output']:
        subdirs.extend(['fast5'])

    if config['dump_adapter_signals']:
        subdirs.extend(['adapter-dumps'])

    for subdir in subdirs:
        fullpath = os.path.join(outputdir, subdir)
        if not os.path.isdir(fullpath):
            os.makedirs(fullpath)

def setup_output_name_mapping(config):
    names = {'fail': OUTPUT_NAME_FAILED}

    if config['barcoding']:
        num_barcodes = config['demultiplexing']['number_of_barcodes']
        for i in range(self.num_barcodes):
            names[i] = OUTPUT_NAME_BARCODES.format(i + 1)
    else:
        names['pass'] = OUTPUT_NAME_PASSED

    if config['filter_unsplit_reads']:
        names['artifact'] = OUTPUT_NAME_ARTIFACT

    return names

def show_configuration(config, args, file):
    _ = partial(print, sep='\t', file=file)
    bool2yn = lambda b: 'Yes' if b else 'No'

    _("== Analysis settings ======================================")
    _(" * Input:", config['inputdir'])
    _(" * Output:", config['outputdir'])
    _(" * Processes:", args.parallel)
    _(" * Presets:", config['preset_name'])
    _(" * Trim 3' adapter:\t", bool2yn(config['trim_adapter']))
    _(" * Filter concatenated read:", bool2yn(config['filter_unsplit_reads']))
    _(" * Separate by barcode:\t", bool2yn(config['barcoding']))
    _(" * Fast5 in output:\t", bool2yn(config['fast5_output']),
           '(Symlink)' if config['fast5_always_symlink'] else '')
    _(" * Basecall table in output:", bool2yn(config['dump_basecalls']))

    if config['dump_adapter_signals']:
        _(" * Dump adapter signals for training:", "Yes")
    _("===========================================================\n")

def main(args):
    if not args.quiet:
        show_banner()

    if not os.path.isdir(args.input):
        errx('ERROR: Cannot open the input directory {}.'.format(args.input))

    if not os.path.isdir(args.output):
        try:
            os.makedirs(args.output)
        except:
            errx('ERROR: Failed to create the output directory {}.'.format(args.output))

    config = load_config(args)
    config['quiet'] = args.quiet
    config['interactive'] = not args.yes
    config['inputdir'] = args.input
    config['outputdir'] = args.output
    config['barcoding'] = args.barcoding
    config['filter_unsplit_reads'] = not args.keep_unsplit
    config['batch_chunk_size'] = args.batch_chunk
    config['dump_adapter_signals'] = args.dump_adapter_signals
    config['dump_basecalls'] = args.dump_basecalled_events
    config['fast5_output'] = args.fast5
    config['fast5_always_symlink'] = args.symlink_fast5
    config['trim_adapter'] = args.trim_adapter
    config['output_names'] = setup_output_name_mapping(config)

    create_output_directories(args.output, config)

    with open(os.path.join(args.output, 'octopus.log'), 'w') as cfgf:
        print("Octopus {}".format(__version__), file=cfgf)
        print("\nStarted at", time.asctime(), file=cfgf)
        print("\nCommand line:", ' '.join(sys.argv) + '\n', file=cfgf)

        show_configuration(config, args, file=cfgf)

    if not config['quiet']:
        show_configuration(config, args, sys.stdout)

    ProcessingSession.run(config, args)


def __main__():
    parser = argparse.ArgumentParser(
        prog='octopus',
        description='A versatile tool for processing nanopore direct RNA sequencing data')

    parser.add_argument('-i', '--input', required=True,
                        help='Path to the directory with the input FAST5 files.')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory path')
    parser.add_argument('-p', '--parallel', default=1, type=int,
                        help='Number of worker processes (default: 1)')
    parser.add_argument('--batch-chunk', default=128, type=int,
                        help='Number of files in a single batch (default: 128)')
    parser.add_argument('-c', '--config', default='',
                        help='Path to signal processing configuration.')
    parser.add_argument('--albacore', default=False, action='store_true',
                        help='Call the ONT albacore for basecalling on-the-fly.')
    parser.add_argument('--barcoding', default=False, action='store_true',
                        help='Sort barcoded reads into separate outputs.')
    parser.add_argument('--trim-adapter', default=False, action='store_true',
                        help="Trim 3' adapter sequences from FASTQ outputs.")
    parser.add_argument('--keep-unsplit', default=False, action='store_true',
                        help="Don't remove unsplit reads fused of two or more RNAs in output.")
    parser.add_argument('--dump-adapter-signals', default=False, action='store_true',
                        help='Dump adapter signal dumps for training')
    parser.add_argument('--dump-basecalled-events', default=False, action='store_true',
                        help='Dump basecalled events from albacore to the output')
    parser.add_argument('-5', '--fast5', default=False, action='store_true',
                        help='Link or copy FAST5 files to separate output directories.')
    parser.add_argument('--symlink-fast5', default=False, action='store_true',
                        help='Create symbolic links to FAST5 files in output directories '
                             'even when hard linking is possible.')
    parser.add_argument('-q', '--quiet', default=False, action='store_true',
                        help='Suppress non-error messages.')
    parser.add_argument('-y', '--yes', default=False, action='store_true',
                        help='Suppress all questions.')

    args = parser.parse_args(sys.argv[1:])
    main(args)
