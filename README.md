# Octopus
Signal-level preprocessor for Oxford Nanopore direct RNA sequencing
(DRS) data. Octopus transforms the data more accessible to RNA Biology.

## Features
* Demultiplexing barcoded direct RNA sequencing libraries
* 3′ Adapter sequence trimming
* Pseudo-fusion read filtering
* Basecalling with ONT `albacore` (even faster than `albacore` itself)
* Live basecalling and processing
* Real-time sequenced read alignments with `minimap2`
* Organizing data in <100 files

## Installation
Octopus requires Python 3.5+ and [pip](http://pypi.python.org/pypi/pip) to install.
This `pip` command installs `octopus` with its essential dependencies. Currently,
you need to install `pomegranate` manually before installing octopus due to the
memory leakage in the released versions of `pomegranate`.

```bash
pip install cython && pip install git+https://github.com/jmschrei/pomegranate.git
pip install octopus
```

To install it together with all optional dependencies (except `albacore`), use this
command:

```bash
pip install 'octopus[full]'
```

### Additional (Optional) Dependency
Octopus requires the FAST5 files basecalled using
[ONT albacore](https://community.nanoporetech.com/downloads) as inputs.
Alternatively, octopus can also internally call `albacore` during the
processing without the prior basecalling if the `albacore` package is
available from the environment.

## Quick Start
Produce FASTQ files with 3′ adapter sequences removed from a directory
containing a bunch of FAST5 files.

```bash
octopus -i path/to/fast5 -o path/to/output --trim-adapter
```

Four barcoded direct RNA sequencing libraries (see below for details)
were pooled and sequenced together. Trim 3′ adapters and demultiplex
the reads into separate FASTQ files.

```bash
octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding
```

In addition to the above, create directories containing hard-links to
the original FAST5 files organized separately by the barcodes.

```bash
octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5
```

In case the FAST5 files are not basecalled yet, just a switch lets
`octopus` call `albacore` internally. Multicore machines help.

```bash
octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5 --albacore-onthefly --parallel 40
```

With the `--live` switch, All tasks can be processed as soon as reads
are produced from MinKNOW.

```bash
octopus -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --albacore-onthefly --parallel 40 --live
```

By default, `octopus` discards pseudo-fusion reads which may originate
from insufficiently segmented signals. You can suppress the filtering
by an option.

```bash
octopus -i path/to/fast5 -o path/to/output --keep-unsplit
```

## Barcoding direct RNA sequencing libraries
ONT direct RNA sequencing libraries are prepared by subsequently attaching
two different 3' adapters, [RTA and RMX](https://community.nanoporetech.com/protocols/sequence-specific-direct-rna-sequencing/v/drss_9035_v1_revg_11may2017/overview-of-the-direct-rna),
respectively. Both are double-stranded DNAs with Y-burged ends on the
3'-sides. Octopus barcoded libraries can be built with modified versions of
RTA adapters. Unlike in the DNA sequencing libraries, `octopus` demultiplexes
in signal-level to ensure the highest accuracy. The distribution contains
demultiplexer models pre-trained with four different DNA barcodes.
Order these sequences and replace the original RTA adapters as many as you
need in the experiment.

```yaml
BC1 Oligo A: 5'-/5Phos/CCTCCCCTAAAAACGAGCCGCATTTGCGTAGTAGGTTC-3'
BC1 Oligo B: 5'-GAGGCGAGCGGTCAATTTTCGCAAATGCGGCTCGTTTTTAGGGGAGGTTTTTTTTTT-3'
```

```yaml
BC2 Oligo A: 5'-/5Phos/CCTCGTCGGTTCTAGGCATCGCGTATGCTAGTAGGTTC-3'
BC2 Oligo B: 5'-GAGGCGAGCGGTCAATTTTGCATACGCGATGCCTAGAACCGACGAGGTTTTTTTTTT-3'
```

```yaml
BC3 Oligo A: 5'-/5Phos/CCTCCCACTTTCACACGCACTAACCAGGTAGTAGGTTC-3'
BC3 Oligo B: 5'-GAGGCGAGCGGTCAATTTTCCTGGTTAGTGCGTGTGAAAGTGGGAGGTTTTTTTTTT-3'
```

```yaml
BC4 Oligo A: 5'-/5Phos/CCTCCTTCAGAAGAGGGTCGCTTCTACCTAGTAGGTTC-3'
BC4 Oligo B: 5'-GAGGCGAGCGGTCAATTTTGGTAGAAGCGACCCTCTTCTGAAGGAGGTTTTTTTTTT-3'
```

## Basecalling with ONT Albacore
Lorem ipsum dolor sit amet.

## Live Basecalling and Processing
Lorem ipsum dolor sit amet.

## Citing Octopus
A pre-print is going to be released soon.
