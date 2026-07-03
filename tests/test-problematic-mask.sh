#!/bin/sh
set -eu

tmp="${TMPDIR:-/tmp}/minibwa-problematic-mask.$$"
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
(tmp / "reads.fa").write_text(">read_masked\n" + read + "\n", encoding="ascii")
(tmp / "mask.bed").write_text("chrMock\t100\t300\tmock_false_dup\n", encoding="ascii")
PY

./minibwa index "$tmp/ref.fa" "$tmp/ref" >/dev/null 2>"$tmp/index.err"
./minibwa map -f "$tmp/ref" "$tmp/reads.fa" >"$tmp/no-mask.paf" 2>"$tmp/no-mask.err"
./minibwa map -f --problematic-bed="$tmp/mask.bed" --problematic-mapq-cap=7 "$tmp/ref" "$tmp/reads.fa" >"$tmp/masked.paf" 2>"$tmp/masked.err"

no_mask_mapq=$(awk '$1 == "read_masked" { print $12; exit }' "$tmp/no-mask.paf")
masked_mapq=$(awk '$1 == "read_masked" { print $12; exit }' "$tmp/masked.paf")

if [ -z "$no_mask_mapq" ] || [ -z "$masked_mapq" ]; then
	echo "missing alignment in problematic-mask fixture" >&2
	exit 1
fi

if [ "$no_mask_mapq" -le 7 ]; then
	echo "expected uncapped MAPQ > 7, got $no_mask_mapq" >&2
	exit 1
fi

if [ "$masked_mapq" -ne 7 ]; then
	echo "expected masked MAPQ cap of 7, got $masked_mapq" >&2
	exit 1
fi

if ! grep -q 'gm:Z:GRC' "$tmp/masked.paf"; then
	echo "expected gm:Z:GRC tag on masked alignment" >&2
	exit 1
fi
