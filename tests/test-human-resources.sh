#!/bin/sh
set -eu

tmp="${TMPDIR:-/tmp}/minibwa-human-resources.$$"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

# Build a small reference plus a GRCh38-named HMF resource bundle.
python3 - "$tmp" <<'PY'
import gzip, pathlib, random, sys

tmp = pathlib.Path(sys.argv[1])
rng = random.Random(7)
def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

(tmp / "ref.fa").write_text(">chr1\n" + seq(2000) + "\n>chr2\n" + seq(2000) + "\n", encoding="ascii")

res = tmp / "bundle"
(res / "dna_pipeline/variants").mkdir(parents=True)
(res / "dna_pipeline/common").mkdir(parents=True)
(res / "dna_pipeline/sv").mkdir(parents=True)
(res / "other/lilac").mkdir(parents=True)

# mappability: 5 intervals of length 50 -> records=5 bases=250 (nested .bed.gz)
with gzip.open(res / "dna_pipeline/variants/mappability_150.38.bed.gz", "wt") as f:
    f.write("".join(f"chr1\t{i*100}\t{i*100+50}\n" for i in range(1, 6)))
# unmap_regions: 2 data rows (50 + 100), one comment header -> records=2 bases=150
(res / "dna_pipeline/common/unmap_regions.38.tsv").write_text(
    "#chrom\tstart\tend\nchr1\t10\t60\nchr2\t20\t120\n", encoding="ascii")
# sv_prep_blacklist: rows 30 + 40 + 10 -> records=3 bases=80
(res / "dna_pipeline/sv/sv_prep_blacklist.38.bed").write_text(
    "chr1\t0\t30\nchr2\t5\t45\nchr2\t50\t60\n", encoding="ascii")
# hla: 1 row 500 -> records=1 bases=500
(res / "other/lilac/hla.38.bed").write_text("chr6\t1000\t1500\n", encoding="ascii")
# germline blacklist: rows 50 + 5 -> records=2 bases=55
(res / "dna_pipeline/variants/KnownBlacklist.germline.38.bed").write_text(
    "chr1\t200\t250\nchr1\t300\t305\n", encoding="ascii")
PY

expect_row() { # id records bases
	line=$(awk -F'\t' -v id="$1" '$1 == id { print $3, $4 }' "$sidecar")
	if [ "$line" != "$2 $3" ]; then
		echo "resource $1: expected records/bases '$2 $3', got '$line'" >&2
		exit 1
	fi
}

check_sidecar() { # sidecar-file
	sidecar="$1"
	if ! grep -q '^#resource_build	GRCh38$' "$sidecar"; then
		echo "expected #resource_build GRCh38 in $sidecar" >&2
		exit 1
	fi
	expect_row mappability_150 5 250
	expect_row unmap_regions 2 150
	expect_row sv_prep_blacklist 3 80
	expect_row hla 1 500
	expect_row known_blacklist_germline 2 55
}

# --- directory import ---
./minibwa index "$tmp/ref.fa" "$tmp/idx_dir" --human-resources "$tmp/bundle" >/dev/null 2>"$tmp/dir.err"
check_sidecar "$tmp/idx_dir.hmf.tsv"

# --- tarball import produces identical stats (incl. nested .bed.gz) ---
( cd "$tmp/bundle" && tar -czf "$tmp/hmftools.tar.gz" . )
./minibwa index "$tmp/ref.fa" "$tmp/idx_tar" --human-resources "$tmp/hmftools.tar.gz" >/dev/null 2>"$tmp/tar.err"
check_sidecar "$tmp/idx_tar.hmf.tsv"
awk -F'\t' 'NR>6 {print $1, $3, $4, $5}' "$tmp/idx_dir.hmf.tsv" > "$tmp/dir.stats"
awk -F'\t' 'NR>6 {print $1, $3, $4, $5}' "$tmp/idx_tar.hmf.tsv" > "$tmp/tar.stats"
if ! diff "$tmp/dir.stats" "$tmp/tar.stats" >/dev/null; then
	echo "tar and directory imports disagree on records/bases/checksums" >&2
	exit 1
fi

# --- missing resources must fail and name what is missing ---
mkdir -p "$tmp/partial"
cp "$tmp/bundle/other/lilac/hla.38.bed" "$tmp/partial/"
if ./minibwa index "$tmp/ref.fa" "$tmp/idx_part" --human-resources "$tmp/partial" >/dev/null 2>"$tmp/part.err"; then
	echo "expected failure on incomplete HMF bundle" >&2
	exit 1
fi
grep -q 'found 1/5 HMF resources' "$tmp/part.err" || { echo "missing found-count error" >&2; exit 1; }
grep -q 'missing mappability_150.38.bed.gz' "$tmp/part.err" || { echo "missing per-resource error" >&2; exit 1; }

# --- mixed GRCh37/GRCh38 builds must be rejected ---
mkdir -p "$tmp/mix"
cp "$tmp/bundle/dna_pipeline/common/unmap_regions.38.tsv" "$tmp/mix/"
cp "$tmp/bundle/dna_pipeline/sv/sv_prep_blacklist.38.bed" "$tmp/mix/"
cp "$tmp/bundle/other/lilac/hla.38.bed" "$tmp/mix/"
cp "$tmp/bundle/dna_pipeline/variants/KnownBlacklist.germline.38.bed" "$tmp/mix/"
python3 - "$tmp" <<'PY'
import gzip, pathlib, sys
tmp = pathlib.Path(sys.argv[1])
with gzip.open(tmp / "mix/mappability_150.37.bed.gz", "wt") as f:
    f.write("chr1\t0\t50\n")
PY
if ./minibwa index "$tmp/ref.fa" "$tmp/idx_mix" --human-resources "$tmp/mix" >/dev/null 2>"$tmp/mix.err"; then
	echo "expected failure on mixed genome builds" >&2
	exit 1
fi
grep -q 'mix genome builds' "$tmp/mix.err" || { echo "missing mixed-build error" >&2; exit 1; }
if [ -e "$tmp/idx_mix.hmf.tsv" ]; then
	echo "sidecar must not be written on mixed-build failure" >&2
	exit 1
fi

# --- duplicate resource in a tar is counted once (no double-count) ---
mkdir -p "$tmp/dup/first" "$tmp/dup/second"
cp "$tmp/bundle/other/lilac/hla.38.bed" "$tmp/dup/first/"
cp "$tmp/bundle/other/lilac/hla.38.bed" "$tmp/dup/second/"
cp "$tmp/bundle/dna_pipeline/common/unmap_regions.38.tsv" "$tmp/dup/"
cp "$tmp/bundle/dna_pipeline/sv/sv_prep_blacklist.38.bed" "$tmp/dup/"
cp "$tmp/bundle/dna_pipeline/variants/KnownBlacklist.germline.38.bed" "$tmp/dup/"
cp "$tmp/bundle/dna_pipeline/variants/mappability_150.38.bed.gz" "$tmp/dup/"
( cd "$tmp/dup" && tar -czf "$tmp/dup.tar.gz" . )
./minibwa index "$tmp/ref.fa" "$tmp/idx_dup" --human-resources "$tmp/dup.tar.gz" >/dev/null 2>"$tmp/dup.err"
check_sidecar "$tmp/idx_dup.hmf.tsv"

echo "human-resources import tests passed"
