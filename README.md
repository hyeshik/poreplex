# Poreplex
Signal-level preprocessor for Oxford Nanopore direct RNA sequencing (DRS) data. 
Poreplex transforms the data into a more accessible format for RNA Biology.

## Features
* Demultiplexing barcoded direct RNA sequencing libraries
* Trimming 3′ adapter sequences
* Filtering pseudo-fusion reads
* Basecalling with ONT `albacore` (even faster than `albacore` itself)
* Live basecalling and processing *during* the sequencing
* Real-time read alignments with `minimap2`
* Full-screen dashboard view for real-time reports
* Organizing data in <100 files

<p align="center">
<img src="https://cowork.narrykim.org/nanopore/octopus/raw/master/doc/images/dashboard.gif" width="480px">
</p>

## Installation
Poreplex requires Python 3.5+ and [pip](http://pypi.python.org/pypi/pip) to install.
This `pip` command installs `poreplex` with its essential dependencies. Currently,
you need to install `pomegranate` manually before installing poreplex due to a
memory leakage issue in the released versions of `pomegranate`. You may use the
following command.

```bash
pip install cython && pip install git+https://github.com/jmschrei/pomegranate.git
pip install poreplex
```

To install it together with all optional dependencies (except `albacore`), use this
command:

```bash
pip install 'poreplex[full]'
```

### Additional (Optional) Dependency
As its inputs, Poreplex requires the FAST5 files that were basecalled using
[ONT albacore](https://community.nanoporetech.com/downloads) in advance.
Alternatively, poreplex can also internally call `albacore` during the
processing without a prior basecalling if the `albacore` package is
available from the environment.

## Quick Start
Produce FASTQ files without 3′ adapter sequences from a bunch of FAST5 files.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter
```

Four direct RNA sequencing libraries can be barcoded, pooled and sequenced 
together. (see below for details) Porplex can demultiplex the librariess into 
separate FASTQ files.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding
```

In addition, Poreplex can create directories containing hard-links to
the original FAST5 files, organized separately by the barcodes.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5
```

In case the FAST5 files are not basecalled yet, just a switch lets
`poreplex` call `albacore` internally. Multicore machines help.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --fast5 --basecall --parallel 40
```

With the `--live` switch, All tasks can be processed as soon as reads
are produced from MinKNOW.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --basecall --parallel 40 --live
```

One may want to output *aligned* reads directly to BAM files instead of
FASTQ outputs. Poreplex streams the processed reads to `minimap2` and update
the BAM outputs real-time. A pre-built index (not a FASTA) generated using
`minimap2` must be provided for this.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --basecall \
  --parallel 40 --live --align GRCz11-transcriptome.mmidx
```

More vibrant feedback is provided if you turn on the dashboard switch.

```bash
poreplex -i path/to/fast5 -o path/to/output --trim-adapter --barcoding --basecall \
  --parallel 40 --live --align GRCz11-transcriptome.mmidx --dashboard
```

By default, `poreplex` discards pseudo-fusion reads which may originate
from insufficiently segmented signals. You can suppress the filtering
by a switch.

```bash
poreplex -i path/to/fast5 -o path/to/output --keep-unsplit
```

## Barcoding direct RNA sequencing libraries
ONT direct RNA sequencing libraries are prepared by subsequently attaching
two different 3' adapters, [RTA and RMX](https://community.nanoporetech.com/protocols/sequence-specific-direct-rna-sequencing/v/drss_9035_v1_revg_11may2017/overview-of-the-direct-rna),
respectively. Both are double-stranded DNAs with Y-burged ends on the
3'-sides. Barcoded libraries for Poreplex can be built with modified versions of
RTA adapters. Unlike in the DNA sequencing libraries, `poreplex` demultiplexes
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

## Basecalling with the ONT Albacore
Most studies requiring signal-level analysis need re-basecalling with
the ONT `albacore` to get the event tables in the FAST5 files. `Poreplex`
can internally call the basecaller core routines of `albacore` to yield
the sequences and tables for the downstream analyses. In fact, running
`albacore` via `poreplex` is remarkably faster than running `albacore`
itself in a multi-core machine thanks to more efficient scheduling of
the computational loads.

<p align="center">
<img src="https://cowork.narrykim.org/nanopore/octopus/raw/master/doc/images/poreplex-albacore-benchmark.jpg" width="520px">
</p>

## Live Basecalling and Processing
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim
veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea
commodo consequat. Duis aute irure dolor in reprehenderit in voluptate
velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat
cupidatat non proident, sunt in culpa qui officia deserunt mollit anim
id est laborum.

## Real-time sequence alignments 
Vitae congue eu consequat ac felis donec et odio pellentesque. Eget
velit aliquet sagittis id consectetur purus. Iaculis at erat
pellentesque adipiscing commodo elit. Nisi scelerisque eu ultrices
vitae auctor. Pellentesque diam volutpat commodo sed egestas egestas
fringilla. Enim nulla aliquet porttitor lacus. Dignissim enim sit amet
venenatis urna. Eros in cursus turpis massa tincidunt dui ut. Quisque
sagittis purus sit amet volutpat. Vitae auctor eu augue ut lectus arcu
bibendum. Eu lobortis elementum nibh tellus molestie nunc non. Sit
amet justo donec enim diam vulputate ut pharetra. Vestibulum lectus
mauris ultrices eros in. Varius quam quisque id diam vel quam elementum
pulvinar etiam. Eget est lorem ipsum dolor sit. Sociis natoque penatibus
et magnis. Quis ipsum suspendisse ultrices gravida dictum fusce. Tellus
integer feugiat scelerisque varius morbi. Egestas quis ipsum suspendisse
ultrices gravida dictum fusce. Nibh tortor id aliquet lectus proin nibh
nisl condimentum id.

## Real-time reports
Interdum posuere lorem ipsum dolor sit. Consectetur a erat nam at lectus
urna duis convallis convallis. Augue eget arcu dictum varius duis at
consectetur lorem donec. Elementum nisi quis eleifend quam adipiscing
vitae proin sagittis nisl. Donec ac odio tempor orci dapibus. Convallis
tellus id interdum velit laoreet. Tempus imperdiet nulla malesuada
pellentesque elit eget. Enim nulla aliquet porttitor lacus luctus. Vel
risus commodo viverra maecenas. Vulputate ut pharetra sit amet aliquam
id. Velit scelerisque in dictum non consectetur a. Id venenatis a
condimentum vitae sapien pellentesque habitant morbi. At tempor commodo
ullamcorper a lacus vestibulum. Lacus vestibulum sed arcu non. Potenti
nullam ac tortor vitae purus faucibus ornare suspendisse sed. Suspendisse
potenti nullam ac tortor vitae purus faucibus ornare. Et leo duis ut
diam quam nulla porttitor massa. Suspendisse faucibus interdum posuere
lorem ipsum. Ut sem viverra aliquet eget sit amet tellus. Semper risus
in hendrerit gravida rutrum quisque non tellus orci.

## Pseudo-fusion filter
Eu turpis egestas pretium aenean pharetra magna ac placerat. Sem
nulla pharetra diam sit amet. Posuere lorem ipsum dolor sit amet
consectetur. Nisl condimentum id venenatis a condimentum vitae
sapien. Eget duis at tellus at urna condimentum mattis pellentesque
id. Odio tempor orci dapibus ultrices in iaculis nunc sed. Leo a diam
sollicitudin tempor. Eget nunc scelerisque viverra mauris. Curabitur vitae
nunc sed velit. Volutpat diam ut venenatis tellus in metus. Nec feugiat
in fermentum posuere urna nec tincidunt praesent semper. Ipsum nunc
aliquet bibendum enim facilisis. Morbi tempus iaculis urna id. Tellus
in metus vulputate eu scelerisque. Eu turpis egestas pretium aenean
pharetra magna ac placerat.

## Output formats

### FASTQ
Mauris pharetra et ultrices neque ornare aenean euismod elementum. Quis
ipsum suspendisse ultrices gravida dictum fusce ut. Pharetra pharetra
massa massa ultricies mi quis hendrerit dolor.

### FAST5
Vel quam elementum pulvinar etiam non quam lacus suspendisse. Sed risus
ultricies tristique nulla. Non sodales neque sodales ut etiam sit amet
nisl purus.

### BAM
Sollicitudin ac orci phasellus egestas. Diam in arcu cursus euismod
quis. Sit amet commodo nulla facilisi. Sollicitudin aliquam ultrices
sagittis orci a scelerisque purus semper.

### Nanopolish database
Leo a diam sollicitudin tempor id eu nisl nunc mi. Euismod quis viverra
nibh cras pulvinar mattis nunc sed. Gravida cum sociis natoque penatibus
et. Ligula ullamcorper malesuada proin libero nunc consequat interdum
varius sit. A iaculis at erat pellentesque adipiscing commodo.

### Basecalled events tables
Quisque id diam vel quam elementum pulvinar etiam non quam. Laoreet
suspendisse interdum consectetur libero id faucibus nisl tincidunt
eget. In est ante in nibh mauris. Pellentesque elit eget gravida cum
sociis natoque penatibus.

## Citing Poreplex
A pre-print is going to be released soon.

