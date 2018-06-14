# Octopus
Octopus is a squiggle-level preprocessor for nanopore direct RNA
sequencing (DRS) data. It transforms data from a sequencing into a more
friendly form to RNA Biology.

## Functions
* Demultiplexes reads from barcoded libraries
* Trims 3′ adapter sequences
* Filters out pseudo-fusion artifacts produced by insufficient gaps
  between two molecules
* Basecalls through the ONT albacore (faster than the original pipeline)
* Provides live basecalling and processing
* Organizes the basecalled data in a less redundant way

## Installation
Octopus requires Python 3.5+ and pip to install.

    pip install octopus

[ONT albacore](https://community.nanoporetech.com/downloads) can be
installed anytime if you want to process the raw FAST5 files with octopus.

## Quick Start
Produce FASTQ files with 3′ adapter sequences removed from a directory
containing a bunch of FAST5 files.

    octopus -i path/to/fast5 -o path/to/output --trim-adapter

Four barcoded direct RNA sequencing libraries (see below for details)
were pooled and sequenced together. Trim 3′ adapters and demultiplex
the reads into separate FASTQ files.

    octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding

In addition to the above, create directories containing hard-links to
the original FAST5 files organized separately by the barcodes.

    octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5

In case the FAST5 files are not basecalled yet, just a switch lets
``octopus`` call albacore internally. Multicore machines help.

    octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5 --albacore-onthefly --parallel 40

By default, ``octopus`` discards pseudo-fusion reads which may originate
from insufficiently segmented signals. You can suppress the filtering
by an option.

    octopus -i path/to/fast5 -o path/to/output --keep-unsplit

## Barcoding direct RNA sequencing libraries
DNA sequences and some explanations.

## Citation
A pre-print is going to released soon.
