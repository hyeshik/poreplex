#
# Copyright (c) 2018-2019 Institute for Basic Science
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
import subprocess as sp
import logging
from functools import partial
from . import *
from .pipeline import ProcessingSession
from .alignment_writer import check_minimap2_index
from .utils import *

VERSION_STRING = """\
poreplex version {version}
Written by Hyeshik Chang <hyeshik@snu.ac.kr>.

Copyright (c) 2018-2019 Institute for Basic Science""".format(version=__version__)

def show_banner():
    print("""
\x1b[1mPoreplex\x1b[0m version {version} by Hyeshik Chang <hyeshik@snu.ac.kr>
- Cuts nanopore direct RNA sequencing data into bite-size pieces for RNA Biology
""".format(version=__version__))

class VersionAction(argparse.Action):

    def __init__(self, option_strings, dest, default=None, required=False,
                 help=None, metavar=None):
        super(VersionAction, self).__init__(
            option_strings=option_strings, dest=dest, nargs=0, help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print(VERSION_STRING)
        parser.exit()

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

def init_logging(config):
    logfile = os.path.join(config['outputdir'], 'poreplex.log')
    logger = logging.getLogger('poreplex')
    handler = logging.FileHandler(logfile, 'w')

    logger.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter('%(asctime)-15s %(message)s'))
    logger.addHandler(handler)

    return logger

def create_output_directories(config):
    outputdir = config['outputdir']
    existing = os.listdir(outputdir)
    if existing:
        while config['interactive']:
            try:
                answer = input('Output directory {} is not empty. Clear it? (y/N) '
                                .format(outputdir))
            except KeyboardInterrupt:
                raise SystemExit
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

    subdirs = []
    conditional_subdirs = [
        ('fastq_output', 'fastq'),
        ('fast5_output', 'fast5'),
        ('nanopolish_output', 'nanopolish'),
        ('minimap2_index', 'bam'),
        ('dump_adapter_signals', 'adapter-dumps'),
        ('dump_basecalls', 'events'),
    ]
    for condition, subdir in conditional_subdirs:
        if config[condition]:
            subdirs.append(subdir)

    for subdir in subdirs:
        fullpath = os.path.join(outputdir, subdir)
        if not os.path.isdir(fullpath):
            os.makedirs(fullpath)

    if not os.path.isdir(config['tmpdir']):
        os.makedirs(config['tmpdir'])
        config['cleanup_tmpdir'] = True


def setup_output_name_mapping(config):
    label_names = {'fail': OUTPUT_NAME_FAILED, 'pass': OUTPUT_NAME_PASSED}

    if config['filter_unsplit_reads']:
        label_names['artifact'] = OUTPUT_NAME_ARTIFACT

    if config['barcoding']:
        num_barcodes = config['demultiplexing']['number_of_barcodes']
        barcode_names = {None: OUTPUT_NAME_UNDETERMINED}
        for i in range(num_barcodes):
            barcode_names[i] = OUTPUT_NAME_BARCODES.format(n=i + 1)

        layout_maps = {
            (label, bc): os.path.join(labelname, bcname)
            for label, labelname in label_names.items()
            for bc, bcname in barcode_names.items()
        }
    else:
        barcode_names = {None: OUTPUT_NAME_BARCODING_OFF}
        layout_maps = {
            (label, None): labelname for label, labelname in label_names.items()}

    return label_names, barcode_names, layout_maps


def show_configuration(config, output):
    if hasattr(output, 'write'): # file-like object
        _ = partial(print, sep='\t', file=output)
    else: # logger object
        _ = lambda *args: output.info(' '.join(map(str, args)))

    bool2yn = lambda b: 'Yes' if b else 'No'

    _("== Analysis settings ======================================")
    _(" * Input:", config['inputdir'],
      '(live, {} sec delay)'.format(config['analysis_start_delay'])
      if config['live'] else '')
    _(" * Output:", config['outputdir'])
    _(" * Processes:", config['parallel'])
    _(" * Presets:", config['preset_name'])
    _(" * Basecall on-the-fly:\t",
        'Yes (albacore {})'.format(config['albacore_version'])
        if config['albacore_onthefly'] else 'No (use previous analyses)')
    _(" * Trim 3' adapter:\t", bool2yn(config['trim_adapter']))
    _(" * Filter concatenated read:", bool2yn(config['filter_unsplit_reads']))
    _(" * Separate by barcode:\t", bool2yn(config['barcoding']))
    _(" * Real-time alignment:\t", bool2yn(config['minimap2_index']))
    _(" * FASTQ in output:\t", bool2yn(config['fastq_output']))
    _(" * FAST5 in output:\t", bool2yn(config['fast5_output']),
           '(Symlink)' if config['fast5_always_symlink'] else '')
    _(" * Basecall table in output:", bool2yn(config['dump_basecalls']))

    if config['dump_adapter_signals']:
        _(" * Dump adapter signals for training:", "Yes")
    _("===========================================================")
    _("")

def test_prerequisite_compatibility(config):
    from distutils.version import LooseVersion
    from pomegranate import __version__ as pomegranate_version
    if LooseVersion(pomegranate_version) <= LooseVersion('0.9.0'):
        errprint('''
WARNING: You have pomegranate {} installed, which has a known
problem that the memory consumption indefinitely grow. The processing
may stop after processing few thousands of reads due to the out of memory
(OOM) errors. Use this command to install until the new release comes out
with the fix:

  pip install cython
  pip install git+https://github.com/jmschrei/pomegranate.git\n'''.format(pomegranate_version))

def test_optional_features(config):
    if config['albacore_onthefly']:
        config['albacore_configuration'] = os.path.join(
            config['outputdir'], 'albacore-configuration.cfg')

        # Check the availability and version compatibility in a subprocess to
        # avoid potential conflicts between duplicated resources in the C++
        # library memory space when the workers are forked into multiple processes.
        result = sp.check_output([sys.executable, '-m',
            'poreplex.basecall_albacore', config['albacore_configuration'],
            config['flowcell'], config['kit']]).decode().strip()
        if result.startswith('okay'):
            config['albacore_version'] = result.split()[1]
        else:
            errx('ERROR: ' + result)

    if config['barcoding']:
        try:
            from .barcoding import BarcodeDemultiplexer
        except:
            errx("ERROR: Barcoding support (--barcoding) requires keras and tensorflow.")

    if config['live']:
        try:
            from inotify.adapters import InotifyTree
        except:
            errx("ERROR: Live monitoring (--live) requires the inotify module.")

def test_inputs_and_outputs(config):
    if not os.path.isdir(config['inputdir']):
        errx('ERROR: Cannot open the input directory {}.'.format(config['inputdir']))

    if not os.path.isdir(config['outputdir']):
        try:
            os.makedirs(config['outputdir'])
        except:
            errx('ERROR: Failed to create the output directory {}.'.format(config['outputdir']))

    if config['minimap2_index']:
        try:
            check_minimap2_index(config['minimap2_index'])
        except:
            errx('ERROR: Could not load a minimap2 index from {}.'.format(config['minimap2_index']))

def fix_options(config):
    printed_any = False

    if config['dashboard'] and not config['minimap2_index']:
        errprint('WARNING: Dashboard is turned off because it is not informative '
                 'without sequence alignments.')
        config['dashboard'] = False
        printed_any = True

    if printed_any:
        errprint('')

def main(args):
    if not args.quiet:
        show_banner()

    config = load_config(args)
    config['quiet'] = args.quiet
    config['interactive'] = not args.yes
    config['parallel'] = args.parallel
    config['inputdir'] = args.input
    config['outputdir'] = args.output
    config['live'] = args.live
    config['analysis_start_delay'] = args.live_delay if args.live else 0
    config['dashboard'] = args.dashboard
    config['contig_aliases'] = args.contig_aliases
    config['tmpdir'] = args.tmpdir if args.tmpdir else os.path.join(args.output, 'tmp')
    config['cleanup_tmpdir'] = False # will be changed during creation of output dirs
    config['barcoding'] = args.barcoding
    config['measure_polya'] = args.polya
    config['filter_unsplit_reads'] = args.filter_chimera
    config['batch_chunk_size'] = args.batch_size
    config['albacore_onthefly'] = args.basecall
    config['dump_adapter_signals'] = args.dump_adapter_signals
    config['dump_basecalls'] = args.dump_basecalled_events
    config['fastq_output'] = args.align is None or args.fastq
    config['fast5_output'] = args.fast5 or args.symlink_fast5 or args.nanopolish
    config['fast5_always_symlink'] = args.symlink_fast5
    config['nanopolish_output'] = args.nanopolish
    config['trim_adapter'] = args.trim_adapter
    config['minimum_sequence_length'] = args.minimum_length
    config['minimap2_index'] = args.align if args.align else None
    config['label_names'], config['barcode_names'], config['output_layout'] = \
        setup_output_name_mapping(config)

    fix_options(config)

    test_inputs_and_outputs(config)
    create_output_directories(config)

    logger = init_logging(config)
    test_prerequisite_compatibility(config)
    test_optional_features(config)

    logger.info('Starting poreplex version {}'.format(__version__))
    logger.info('Command line: ' + ' '.join(sys.argv))

    show_configuration(config, output=logger)
    if not config['quiet']:
        show_configuration(config, output=sys.stdout)

    procresult = ProcessingSession.run(config, logger)

    if procresult is not None:
        if not config['quiet']:
            procresult(sys.stdout)
        procresult(logger)

    logger.info('Finished.')

    if config['cleanup_tmpdir']:
        try:
            shutil.rmtree(config['tmpdir'])
        except:
            pass

def __main__():
    parser = argparse.ArgumentParser(
        prog='poreplex', add_help=False,
        description='Cuts nanopore direct RNA sequencing data '
                    'into bite-size pieces for RNA Biology')

    group = parser.add_argument_group('Data Settings')
    group.add_argument('-i', '--input', required=True, metavar='DIR',
                       help='path to the directory with the input FAST5 files '
                            '(Required)')
    group.add_argument('-o', '--output', required=True, metavar='DIR',
                       help='output directory path (Required)')
    group.add_argument('-c', '--config', default='', metavar='NAME',
                       help='path to signal processing configuration')

    group = parser.add_argument_group('Basic Processing Options')
    group.add_argument('--trim-adapter', default=False, action='store_true',
                       help="trim 3' adapter sequences from FASTQ outputs")
    group.add_argument('--minimum-length', default=10, type=int, metavar='LEN',
                       help="discard reads shorter than LEN (default: 10)")
    group.add_argument('--filter-chimera', default=False, action='store_true',
                       help="remove unsplit reads fused of two or more RNAs in output")

    group = parser.add_argument_group('Optional Analyses')
    group.add_argument('--barcoding', default=False, action='store_true',
                       help='sort barcoded reads into separate outputs')
    group.add_argument('--polya', default=False, action='store_true',
                       help='output poly(A) tail length measurements')
    group.add_argument('--basecall', default=False, action='store_true',
                       help='call the ONT albacore for basecalling on-the-fly')
    group.add_argument('--align', default=None, type=str, metavar='INDEXFILE',
                       help='align basecalled reads using minimap2 and create BAM files')

    group = parser.add_argument_group('Live Mode')
    group.add_argument('--live', default=False, action='store_true',
                       help='monitor new files in the input directory')
    group.add_argument('--live-delay', default=60, type=int, metavar='SECONDS',
                       help='time to delay the start of analysis in live mode '
                            '(default: 60)')

    group = parser.add_argument_group('Output Options')
    group.add_argument('--fastq', default=False, action='store_true',
                       help='write to FASTQ files even when BAM files are produced')
    group.add_argument('--fast5', default=False, action='store_true',
                       help='link or copy FAST5 files to separate output directories')
    group.add_argument('--symlink-fast5', default=False, action='store_true',
                       help='create symbolic links to FAST5 files in output directories '
                            'even when hard linking is possible')
    group.add_argument('--nanopolish', default=False, action='store_true',
                       help='create a nanopolish readdb to enable access from nanopolish')
    group.add_argument('--dump-adapter-signals', default=False, action='store_true',
                       help='dump adapter signal dumps for training')
    group.add_argument('--dump-basecalled-events', default=False, action='store_true',
                       help='dump basecalled events to the output')

    group = parser.add_argument_group('User Interface')
    group.add_argument('--dashboard', dest='dashboard', default=False, action='store_true',
                       help="show the full screen dashboard")
    group.add_argument('--contig-aliases', default=None, metavar='FILE', type=str,
                       help='path to a tab-separated text file for aliases to show '
                            'as a contig names in the dashboard (see README)')
    group.add_argument('-q', '--quiet', default=False, action='store_true',
                       help='suppress non-error messages')
    group.add_argument('-y', '--yes', default=False, action='store_true',
                       help='suppress all questions')

    group = parser.add_argument_group('Pipeline Options')
    group.add_argument('-p', '--parallel', default=1, type=int, metavar='COUNT',
                       help='number of worker processes (default: 1)')
    group.add_argument('--tmpdir', default='', type=str, metavar='DIR',
                       help='temporary directory for intermediate data')
    group.add_argument('--batch-size', default=128, type=int, metavar='SIZE',
                       help='number of files in a single batch (default: 128)')
    group.add_argument('--version', action=VersionAction,
                       help="show program's version number and exit")
    group.add_argument('-h', '--help', action='help',
                       help='show this help message and exit')

    args = parser.parse_args(sys.argv[1:])
    main(args)
