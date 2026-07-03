[![Build Status](https://github.com/lh3/minibwa/actions/workflows/build.yml/badge.svg)](https://github.com/lh3/minibwa/actions)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/minibwa.svg?style=flag&label=Bioconda)](https://anaconda.org/bioconda/minibwa)
[![Homebrew](https://img.shields.io/homebrew/v/minibwa)](https://formulae.brew.sh/formula/minibwa)
[![preprint](https://img.shields.io/badge/arXiv-2606.15357-blue)](https://arxiv.org/abs/2606.15357)

## Getting Started
```sh
git clone https://github.com/lh3/minibwa
cd minibwa && make

# with test data
./minibwa index tests/data/chrM-human.fa.gz chrM-human              # index the genome
./minibwa map chrM-human tests/data/chrM-read_?.fa.gz > aln.sam     # align and output in SAM

# other examples without test data
minibwa map -ft16 ref.index long-read.fq > aln.paf            # align long reads
minibwa map --hic ref.index reads.interleaved.fq > aln.sam    # align Hi-C short reads

# align *directional* bisulfite sequencing (BS-seq) reads
minibwa index --meth -t8 ref.fa                               # generate BS-seq index
minibwa map --meth ref.fa read1.fq read2.fq > aln.sam         # map BS-seq reads
```

## Introduction

Minibwa aligns short reads against a reference genome. It is the successor of
[bwa-mem][bwa] with a different algorithm. Minibwa is over three times as fast as the
original bwa-mem and twice as fast as [bwa-mem2][bwa-mem2] at comparable accuracy. While
minibwa works with accurate long reads, [minimap2][mm2] is more robust under high
error rate.

Minibwa is a hybrid of bwa-mem and minimap2: it indexes the genome with
Burrow-Wheeler Transform (BWT), finds variable-length seeds like bwa-mem, and
performs chaining and SIMD-based nucleotide alignment with the minimap2
algorithm. Minibwa speeds up bwa-mem2 further with additional prefetch for
seeding, new heuristics to skip unnecessary mate rescue and reduced effort in
highly repetitive regions where reads would often be wrongly mapped due to
structural changes anyway.

## Users' Guide

### Intended use cases

Minibwa is designed for mapping short reads and accurate long reads. It does
not support spliced alignment and has not been tuned for aligning long contigs.
For now, minibwa does not properly work with alternate contigs in the reference
genome. Please use a version of the reference without such contigs.

### Installation

Minibwa requires either SSE4.2 on x86 CPUs or NEON on ARM. It depends on
[zlib][zlib] installed on your system and also includes slightly modified
source code of [mimalloc][mimalloc] and [libsais][libsais] which optionally
uses OpenMP for multi-threading. You can build minibwa with
```sh
make             # automatically detect OpenMP and arm64 vs. x86_64
make omp=0       # disable multi-threading in libsais (no effect on mapping)
make gpl=0       # disable GPL'd code for low-memory BWT building (no effect on mapping)
make mimalloc=0  # disable mimalloc and use the system malloc+kalloc instead
```
This produces a single binary `minibwa` which you can copy to your `PATH`.

For deployment builds on CPU-bound mapping workloads, the Makefile also
provides explicit optimized release paths without changing the default developer
build:
```sh
make release-lto

make pgo-generate
PGO_TRAIN_CMD='./minibwa map -t1 -P -o /dev/null ref.index reads.fastq' make pgo-train
make pgo-use
```
Use a representative indexed reference and read set for `PGO_TRAIN_CMD`; mapping
training data should resemble the production workload you want to optimize.
`make pgo` runs the three PGO stages in order when `PGO_TRAIN_CMD` is supplied.
Clang builds require `llvm-profdata`; `make pgo-use` merges raw profiles before
rebuilding with them.
Advanced builds can also use `make lto=1`, `make pgo=generate`, or
`make pgo=use` directly.

Minibwa also supports CMake for package-manager and IDE integrations:
```sh
cmake -S . -B build/cmake
cmake --build build/cmake
```

### Usage

Like bwa-mem, minibwa requires to index the genome before read alignment.

#### Indexing

You can index the reference genome with
```sh
minibwa index -t8 ref.fa     # index with 8 threads, using 18N RAM (N is the genome size)
minibwa index ref.fa prefix  # use a different index prefix instead of ref.fa
minibwa index -l ref.fa      # use less memory at the cost of performance
minibwa index --meth ref.fa  # generate BS-seq index
minibwa index --human ref.fa # prevent seeds through ambiguous reference bases
minibwa index ref.fa prefix --human-resources hmftools.tar.gz  # add HMF resource manifest
```
Minibwa generates two files: `ref.fa.l2b` for 2-bit encoded reference genome
sequences and `ref.fa.mbw` for BWT and sampled suffix array. In the `--meth`
mode, minibwa additionally generates `ref.fa.meth.mbw` for the BWT of the
3-base genome. With `--human-resources`, minibwa reads the HMF mappability,
unmap-region, blacklist, and HLA BED/TSV resources from an extracted directory
or `hmftools` tarball and writes `ref.fa.hmf.tsv` with resource counts,
checksums, and genome-build metadata.

Lowercase bases in the reference FASTA are treated as soft-masked and their
intervals are recorded in `ref.fa.l2b`. Soft-masking does not change the indexed
sequence, so a reference without lowercase bases behaves exactly as before.

#### Mapping

By default, minibwa dynamically changes multiple internal parameters based on
individual read lengths. It works for both short and accurate long reads.
```sh
minibwa map -t8 ref.fa read1.fq read2.fq    # map paired-end reads and output SAM
minibwa map --single-end -t8 ref.fa read.fq # map R1-only/single-end short reads without pairing work
minibwa map -ft8 ref.fa read.fa.gz          # map single-end or long reads; output PAF
minibwa map --hic ref.fa hic1.fq hic2.fq    # map Hi-C short reads
minibwa map --meth ref.fa read1.fq read2.fq # map BS-seq reads; requiring "index --meth"
minibwa map --human ref.fa read1.fq read2.fq # ALT/HLA-aware SAM XA/SA tags
minibwa map --human-tags ref.fa read1.fq read2.fq # add zc/zm/zh human tags
minibwa map --numt ref.fa read.fq           # cap/tag ambiguous chrM-vs-nuclear hits
minibwa map --human-profile=hmf-grch38 ref.fa read.fq # HMF human policies/resources
```
Use `--single-end` (equivalent to `-P`) for R1-only or other confirmed
single-end short-read inputs when pairing and mate rescue are not needed. This
keeps the same read mapping behavior while skipping paired-end statistics,
pair selection, and mate rescue work.

With `--human` (or `-j` in the `mem` subcommand), the SAM `XA` and `SA` tags
become aware of human ALT/HLA contigs: hits on ALT or HLA contigs that compete
with a primary-assembly hit at the same query locus are reported as alternative
placements in the `XA` tag instead of split-alignment `SA` tags. Exact-score
ties are also resolved toward primary human coordinates before ALT, HLA, decoy,
random, unplaced, and viral/phix contigs.

With `--human-profile=hmf-grch38`, minibwa enables an HMF-compatible HLA policy:
exact score ties between an HLA allele/ALT contig and the main chr6 placement
prefer chr6 as the primary alignment. This avoids primary HLA ALT-contig records
that are incompatible with downstream HMF/oncoanalyser HLA typing and variant
calling. Use `--hla-policy=allele-contig` to preserve allele-contig primary
output; stronger allele-contig alignments are not demoted by the default policy.

With `--human-tags`, minibwa adds optional human annotations to each mapped SAM
and PAF record: `zc:Z` gives the contig class inferred from the contig name
(`primary`, `alt`, `hla`, `decoy`, `random`, `unplaced`, or `unknown`), `zm:f`
gives the unmasked fraction of the alignment computed from soft-masked reference
intervals, and `zh:Z` gives a confidence class (`high_confidence`,
`low_mappability`, `non_primary`, or `unknown`). The SAM header describes these
tags in `@CO` comment lines.

For human references, `minibwa map --mapq-mappability mappability_150.38.bed.gz ref.fa reads.fq`
loads a BED or bedGraph mappability track once and caps primary MAPQ for
alignments overlapping low-mappability intervals. Three-column BED rows use
`--mapq-low-cap` (20 by default); numeric fourth-column scores in the range
0..1 are treated as mappability fractions and converted to a 0..60 cap.

Note in the default adaptive mode, `-g`/`-w`/`-W`/`-N`/`-m`/`-s` only changes
the short-read setting; the long-read setting is fixed. This mode is disabled
with `--adap=no` or when `-x sr` or `-x lr` is specified.

#### Pseudoautosomal regions (chrX/chrY)

When the reference is recognized as a human GRCh37 or GRCh38 assembly, minibwa
detects the built-in chrX/chrY pseudoautosomal regions (PAR1 and PAR2) from the
loaded reference dictionary. Reads placed at equivalent chrX and chrY PAR
coordinates are grouped so the chrX/chrY twin does not collapse the mapping
quality; minibwa keeps a deterministic chrX hit as the primary and tags PAR
alignments with a `pa:Z:PAR1` or `pa:Z:PAR2` tag in the SAM/PAF output.
References that are not recognized as GRCh37 or GRCh38 keep the ordinary
minibwa behavior.

If the reference was indexed with soft-masked (lowercase) intervals, minibwa
uses them to reduce spurious placements: high-occurrence anchors falling
entirely within a soft-masked interval are dropped during seeding, and the
mapping quality of an alignment is capped at 10 when most of its reference span
is soft-masked. References with no soft-masked bases are unaffected.

For GRCh38, `--grch38-mask` caps and tags alignments that overlap the compact
GRC problematic/false-duplication exclusion mask bundled with minibwa. Use
`--problematic-mapq-cap=INT` to change the default cap of 10. You can also
provide a BED file explicitly with `--problematic-bed=FILE`; overlapping SAM/PAF
records get the local tag `gm:Z:GRC` and remain in the output.

#### HMF unmap regions (GRCh38)

Minibwa can annotate hits that fall inside the Hartwig Medical Foundation (HMF)
`unmap_regions` for GRCh38:
```sh
minibwa map --human-profile=hmf-grch38 ref.fa read1.fq read2.fq   # auto-discover unmap_regions.38.tsv
minibwa map --unmap-regions=unmap_regions.38.tsv ref.fa reads.fq  # use an explicit file
```
`--human-profile=hmf-grch38` enables the HMF-compatible HLA primary policy and
looks for the standard `unmap_regions.38.tsv` resource in the working directory,
in the `MINIBWA_HMF_RESOURCES`, `HMF_RESOURCES`, or `HMF_RESOURCE_DIR`
directories, and next to the index prefix. Use `--unmap-regions=FILE` to point
at the TSV explicitly (the file has
`Chromosome`, `PosStart`, `PosEnd`, and `MaxDepth` columns with 1-based
inclusive coordinates). Hits overlapping an unmap region are tagged with
`ur:Z:unmap` and `ud:i:<depth>` (the maximum depth of the overlapping regions)
in the SAM/PAF output. This only adds annotation tags; it does not change read
placement or mapping quality.

For human references with both mitochondrial and nuclear contigs, `--numt`
enables a conservative chrM-vs-nuclear policy. Exact score ties between chrM
and nuclear hits are made deterministic by preferring the nuclear placement,
near-tie chrM/nuclear competitors are capped at MAPQ 10, and ambiguous hits are
annotated with `ng:Z:chrM-nuclear` in SAM or PAF output. Reads unique to chrM
are not affected.

#### Mapping with legacy bwa-mem CLI

Minibwa also provides legacy bwa-mem command-line interface (CLI) via the `mem` subcommand.
However, due to algorithm and parameter differences, many bwa-mem options are ignored.
The output minibwa alignment is also not identical to bwa-mem.

## Developers' Guide

Minibwa provides basic APIs for loading index and aligning reads.
[examples/ex-one.c](examples/ex-one.c) shows an example to align each read
independently; [examples/ex-batch.c](examples/ex-batch.c) aligns multiple reads
in batch, which is faster and also supports paired-end mapping. Run `make
examples` to build both examples. Run `make test` to build minibwa and run the
`tests/` shell tests, including the problematic-region mask and single-end
fast-path regressions.
[docs/dev.md](docs/dev.md) explains how minibwa differs from BWA-MEM and minimap2.

### Project layout

Project-owned implementation lives in [src](src), public headers in
[include](include), bundled dependencies in [third_party](third_party), API
examples in [examples](examples), test data in [tests/data](tests/data), and
paper/manpage sources in [docs](docs).

### Benchmarking

The [bench/](bench/) directory contains helper scripts for measuring mapping
performance. `bench/run-local-bench.sh [out.tsv]` builds minibwa, times indexing,
paired-end mapping, R1-only mapping with and without `--single-end`, and the
`bench` subcommand, then writes a TSV of timings
(default `.context/local-bench.tsv`). It uses the bundled `tests/data` inputs by
default; set `REF_FASTA`, `READ1`, `READ2`, or `BENCH_ITERATIONS` to benchmark
other data.
`bench/run-remote-bench.sh <ssh-host>` copies the source tree to a remote SSH
host, rebuilds there, and optionally times the `bench` subcommand when
`BENCH_INDEX` points at a remote `.mbw` index. `REMOTE_BASE`, `REMOTE_WORK`,
`BENCH_ITERATIONS`, and `LOCAL_LOG` customize paths and run size; logs are
written under `.context/` by default.
`bench/compare-upstream.sh` clones/builds `https://github.com/lh3/minibwa`,
checks that shared paired-end and single-end SAM bodies match after ignoring the
`@PG` command line, and runs repeated timing comparisons. Set `REPEAT` and
`SLOWDOWN_LIMIT` to tune the performance gate, or run the same check with
`make upstream-compare`.
`bench/run-human-benchmark.py` adds a human-focused QA harness. With no inputs it
generates a small GRCh38-shaped synthetic fixture covering primary, HLA-like,
alternate, decoy, and low-mappability placements; `make test` runs this CI-sized
check. The same script can collect metrics for prepared HMF/GIAB runs with
external FASTA, FASTQ, and BED stratification inputs. See
[docs/human-benchmark.md](docs/human-benchmark.md).

## License

Minibwa is distributed under the MIT license. It also incorporates source code
from the following projects:

 * libsais: Apache 2 License. Copyright (c) 2021-2025 Ilya Grebnov
 * mimalloc: MIT License. Copyright (c) 2018-2026 Microsoft Corporation, Daan Leijen

The master branch is optionally built on the following projects:

 * QSufSort: HPND License. Copyright (c) 1999 N. Jesper Larsson
 * bwtgen: GPL 2 License. Copyright (c) 2004 Wong Chi Kwong

Notably, the master branch includes GPL'd [bwtgen.c](third_party/bwtgen/bwtgen.c)
for low-memory BWT construction. If you compile this file, which is the default,
the resulting binary will be GPL'd. You can disable the low-memory algorithm
with `make gpl=0` or `cmake -DMINIBWA_ENABLE_GPL=OFF` to generate non-GPL
binary. The [Apache2 branch][apache2] does not include GPL'd source code.

## Limitations

* Minibwa does not work with noisy long reads or spliced RNA-seq reads.
* Minibwa does not support undirectional bisulfite sequencing data.
* Minibwa does not recognize alternate haplotypes.

[apache2]: https://github.com/lh3/minibwa/tree/Apache2
[zlib]: https://zlib.net/
[mimalloc]: https://github.com/microsoft/mimalloc
[libsais]: https://github.com/IlyaGrebnov/libsais
[bwa]: https://github.com/lh3/bwa
[mm2]: https://github.com/lh3/minimap2
[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
