# Human Benchmark Harness

`bench/run-human-benchmark.py` is a lightweight benchmark and QA harness for the
Human Mapping Workload. It records primary-coordinate stability, MAPQ
distribution, alignment rate, runtime, peak RSS, ALT/HLA/decoy primary counts,
and optional BED-stratified primary counts.

## Synthetic CI Fixture

Run:

```sh
make test
```

or directly:

```sh
./bench/run-human-benchmark.py --out-dir .context/human-benchmark
```

The default run builds a small synthetic reference with:

- a Primary Assembly contig (`chr1`);
- an HLA-like primary contig (`chr6`);
- an Alternate Contig (`chr6_GL000250v2_alt`);
- a Decoy Sequence (`hs37d5_decoy`);
- duplicated low-mappability sequence shared by `chr1` and `chr6`.

The harness writes:

- `.context/human-benchmark/aln.sam`
- `.context/human-benchmark/map.stderr`
- `.context/human-benchmark/metrics.json`
- `.context/human-benchmark/metrics.tsv`

The synthetic run validates stable primary coordinates for unique primary,
alternate, and decoy reads, verifies low MAPQ for a duplicated read, and requires
the synthetic reads to map.

## Prepared HMF/GIAB Runs

On a machine with the HMF GRCh38 bundle staged, run the same harness against
real inputs. Example:

```sh
make
./bench/run-human-benchmark.py \
  --reference /data/hmf/GRCh38.fa.gz \
  --reads1 /data/hmf/smoke/read_1.fq.gz \
  --reads2 /data/hmf/smoke/read_2.fq.gz \
  --bed /data/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  --bed /data/hmf/stratifications/low_mappability.bed \
  --bed /data/hmf/stratifications/hla.bed \
  --bed /data/hmf/stratifications/blacklist.bed \
  --out-dir .context/hmf-giab-human-benchmark \
  --threads 8 \
  --no-validate
```

Use `--index-prefix /path/to/existing/prefix` to reuse an existing `.mbw`/`.l2b`
index. Add `--expect-primary-tsv expected.tsv` when a targeted truth set has
known placements. The TSV format is:

```text
read_name	contig	one_based_pos	max_distance
```

The fourth column is optional and defaults to zero. External runs with
`--no-validate` collect metrics without enforcing thresholds, which is useful
for first-pass smoke tests on `do-sydney` or another prepared host.
