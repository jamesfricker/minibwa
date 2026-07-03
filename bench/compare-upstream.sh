#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-"$ROOT/.context/upstream-compare.tsv"}
UPSTREAM_DIR=${UPSTREAM_DIR:-"$ROOT/.context/upstream-minibwa"}
UPSTREAM_URL=${UPSTREAM_URL:-"https://github.com/lh3/minibwa.git"}
REF_FASTA=${REF_FASTA:-"$ROOT/tests/data/chrM-human.fa.gz"}
READ1=${READ1:-"$ROOT/tests/data/chrM-read_1.fa.gz"}
READ2=${READ2:-"$ROOT/tests/data/chrM-read_2.fa.gz"}
REPEAT=${REPEAT:-7}
SLOWDOWN_LIMIT=${SLOWDOWN_LIMIT:-1.25}

for input in "$REF_FASTA" "$READ1" "$READ2"; do
	if [ ! -f "$input" ]; then
		printf "missing input: %s\n" "$input" >&2
		exit 1
	fi
done

if [ ! -d "$UPSTREAM_DIR/.git" ]; then
	mkdir -p "$(dirname "$UPSTREAM_DIR")"
	git clone --depth 1 "$UPSTREAM_URL" "$UPSTREAM_DIR"
fi

mkdir -p "$ROOT/.context"
TMP=$(mktemp -d "$ROOT/.context/upstream-compare.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

cd "$ROOT"
make >/dev/null
make -C "$UPSTREAM_DIR" >/dev/null

CURRENT="$ROOT/minibwa"
BASE="$UPSTREAM_DIR/minibwa"

"$CURRENT" index "$REF_FASTA" "$TMP/current-ref" >"$TMP/current-index.out" 2>"$TMP/current-index.err"
"$BASE" index "$REF_FASTA" "$TMP/base-ref" >"$TMP/base-index.out" 2>"$TMP/base-index.err"

run_map() {
	label=$1
	bin=$2
	shift 2
	"$bin" "$@" >"$TMP/$label.sam" 2>"$TMP/$label.err"
	grep -v '^@PG' "$TMP/$label.sam" >"$TMP/$label.body"
}

run_map current_pe "$CURRENT" map -t1 -x sr "$TMP/current-ref" "$READ1" "$READ2"
run_map base_pe "$BASE" map -t1 -x sr "$TMP/base-ref" "$READ1" "$READ2"
run_map current_single "$CURRENT" map -t1 -x sr --single-end "$TMP/current-ref" "$READ1"
run_map base_single "$BASE" map -t1 -x sr -P "$TMP/base-ref" "$READ1"

cmp "$TMP/current_pe.body" "$TMP/base_pe.body"
cmp "$TMP/current_single.body" "$TMP/base_single.body"

printf "label\tbase_real_seconds\tcurrent_real_seconds\tcurrent_over_base\n" >"$OUT"

time_command() {
	label=$1
	bin=$2
	shift 2
	i=0
	while [ "$i" -lt "$REPEAT" ]; do
		"$bin" "$@" >/dev/null 2>"$TMP/time-$label-$i.err"
		awk '/Real time:/ {
			for (j = 1; j <= NF; ++j)
				if ($j == "time:") print $(j + 1)
		}' "$TMP/time-$label-$i.err"
		i=$((i + 1))
	done | awk '{ sum += $1; n += 1 } END { if (n == 0) exit 1; printf "%.6f\n", sum / n }'
}

record_perf() {
	label=$1
	base_mean=$2
	current_mean=$3
	ratio=$(awk -v c="$current_mean" -v b="$base_mean" 'BEGIN { printf "%.6f", c / b }')
	printf "%s\t%s\t%s\t%s\n" "$label" "$base_mean" "$current_mean" "$ratio" >>"$OUT"
	awk -v ratio="$ratio" -v limit="$SLOWDOWN_LIMIT" -v label="$label" '
		BEGIN {
			if (ratio > limit) {
				printf "performance regression for %s: current/base %.3f > %.3f\n", label, ratio, limit > "/dev/stderr"
				exit 1
			}
		}
	'
}

base_mean=$(time_command "base-map-paired" "$BASE" map -t1 -x sr "$TMP/base-ref" "$READ1" "$READ2")
current_mean=$(time_command "current-map-paired" "$CURRENT" map -t1 -x sr "$TMP/current-ref" "$READ1" "$READ2")
record_perf map_paired "$base_mean" "$current_mean"

base_mean=$(time_command "base-map-single-fast" "$BASE" map -t1 -x sr -P "$TMP/base-ref" "$READ1")
current_mean=$(time_command "current-map-single-fast" "$CURRENT" map -t1 -x sr --single-end "$TMP/current-ref" "$READ1")
record_perf map_single_fast "$base_mean" "$current_mean"

printf "upstream comparison passed; wrote %s\n" "$OUT"
