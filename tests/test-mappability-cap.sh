#!/bin/sh
set -eu

tmp="${TMPDIR:-/tmp}/minibwa-mappability-cap.$$"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" <<'PY'
import pathlib
import random
import sys

tmp = pathlib.Path(sys.argv[1])
rng = random.Random(23)
seq = "".join(rng.choice("ACGT") for _ in range(600))
read = seq[120:270]

(tmp / "ref.fa").write_text(">chrMock\n" + seq + "\n", encoding="ascii")
(tmp / "reads.fa").write_text(">read_lowmap\n" + read + "\n", encoding="ascii")
# A SINGLE low-mappability interval covering the read's placement. This exercises
# the single-interval sort_and_merge path where max_en must be initialised, or the
# cap is silently never applied.
(tmp / "map.bed").write_text("chrMock\t100\t300\n", encoding="ascii")
# bedGraph with a fractional mappability score of 0.10 -> cap round(60*0.10) = 6.
(tmp / "map.bedgraph").write_text("chrMock\t100\t300\t0.10\n", encoding="ascii")
PY

./minibwa index "$tmp/ref.fa" "$tmp/ref" >/dev/null 2>"$tmp/index.err"
./minibwa map -f "$tmp/ref" "$tmp/reads.fa" >"$tmp/no-cap.paf" 2>"$tmp/no-cap.err"
./minibwa map -f --mapq-mappability "$tmp/map.bed" --mapq-low-cap 5 "$tmp/ref" "$tmp/reads.fa" >"$tmp/capped.paf" 2>"$tmp/capped.err"
./minibwa map -f --mapq-mappability "$tmp/map.bedgraph" "$tmp/ref" "$tmp/reads.fa" >"$tmp/frac.paf" 2>"$tmp/frac.err"

no_cap_mapq=$(awk '$1 == "read_lowmap" { print $12; exit }' "$tmp/no-cap.paf")
capped_mapq=$(awk '$1 == "read_lowmap" { print $12; exit }' "$tmp/capped.paf")
frac_mapq=$(awk '$1 == "read_lowmap" { print $12; exit }' "$tmp/frac.paf")

if [ -z "$no_cap_mapq" ] || [ -z "$capped_mapq" ] || [ -z "$frac_mapq" ]; then
	echo "missing alignment in mappability-cap fixture" >&2
	exit 1
fi

if [ "$no_cap_mapq" -le 5 ]; then
	echo "expected uncapped MAPQ > 5, got $no_cap_mapq" >&2
	exit 1
fi

# Single 3-column interval must apply the --mapq-low-cap value.
if [ "$capped_mapq" -ne 5 ]; then
	echo "expected single-interval low-mappability MAPQ cap of 5, got $capped_mapq" >&2
	exit 1
fi

# bedGraph fractional score 0.10 must convert to a 0..60 cap of 6.
if [ "$frac_mapq" -ne 6 ]; then
	echo "expected fractional mappability MAPQ cap of 6, got $frac_mapq" >&2
	exit 1
fi

echo "mappability-cap test passed: baseline=$no_cap_mapq low-cap=$capped_mapq fractional=$frac_mapq"
