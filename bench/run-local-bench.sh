#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-"$ROOT/.context/local-bench.tsv"}
mkdir -p "$ROOT/.context"
TMP=$(mktemp -d "$ROOT/.context/local-bench.XXXXXX")
DATA_DIR=${DATA_DIR:-"$ROOT/tests/data"}
REF_FASTA=${REF_FASTA:-"$DATA_DIR/chrM-human.fa.gz"}
READ1=${READ1:-"$DATA_DIR/chrM-read_1.fa.gz"}
READ2=${READ2:-"$DATA_DIR/chrM-read_2.fa.gz"}
BENCH_ITERATIONS=${BENCH_ITERATIONS:-100000}
INDEX_PREFIX="$TMP/reference"
INDEX_FILE="$INDEX_PREFIX.mbw"

for input in "$REF_FASTA" "$READ1" "$READ2"; do
	if [ ! -f "$input" ]; then
		printf "missing input: %s\n" "$input" >&2
		exit 1
	fi
done

cd "$ROOT"
make >/dev/null

printf "timestamp\thost\tgit_rev\tlabel\treal_seconds\tcpu_seconds\tpeak_rss_gb\textra\n" >"$OUT"

timestamp() { date -u +"%Y-%m-%dT%H:%M:%SZ"; }
git_rev() { git rev-parse --short HEAD 2>/dev/null || printf "unknown"; }

append_minibwa_timing() {
	label=$1
	log=$2
	extra=$3
	awk -v ts="$(timestamp)" -v host="$(hostname)" -v rev="$(git_rev)" -v label="$label" -v extra="$extra" '
		/Real time:/ {
			for (i = 1; i <= NF; ++i) {
				if ($i == "time:") real = $(i + 1)
				else if ($i == "CPU:") cpu = $(i + 1)
				else if ($i == "RSS:") rss = $(i + 1)
			}
		}
		END { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", ts, host, rev, label, real, cpu, rss, extra }
	' "$log" >>"$OUT"
}

append_bwt_timing() {
	label=$1
	log=$2
	extra=$3
	awk -v ts="$(timestamp)" -v host="$(hostname)" -v rev="$(git_rev)" -v label="$label" -v extra="$extra" '
		/^t = / { cpu = $3 }
		/checksum = / { checksum = $3 }
		END { printf "%s\t%s\t%s\t%s\t\t%s\t\t%s checksum=%s\n", ts, host, rev, label, cpu, extra, checksum }
	' "$log" >>"$OUT"
}

./minibwa index "$REF_FASTA" "$INDEX_PREFIX" >"$TMP/index.out" 2>"$TMP/index.err"
append_minibwa_timing "index_reference" "$TMP/index.err" "ref=$REF_FASTA"

./minibwa map "$INDEX_PREFIX" "$READ1" "$READ2" >"$TMP/aln.sam" 2>"$TMP/map.err"
append_minibwa_timing "map_paired_sam" "$TMP/map.err" "reads=$READ1,$READ2 sam_bytes=$(wc -c < "$TMP/aln.sam" | tr -d ' ')"

./minibwa map "$INDEX_PREFIX" "$READ1" >"$TMP/aln-single-default.sam" 2>"$TMP/map-single-default.err"
append_minibwa_timing "map_single_sam_default" "$TMP/map-single-default.err" "reads=$READ1 sam_bytes=$(wc -c < "$TMP/aln-single-default.sam" | tr -d ' ')"

./minibwa map --single-end "$INDEX_PREFIX" "$READ1" >"$TMP/aln-single-fast.sam" 2>"$TMP/map-single-fast.err"
append_minibwa_timing "map_single_sam_fast" "$TMP/map-single-fast.err" "reads=$READ1 option=--single-end sam_bytes=$(wc -c < "$TMP/aln-single-fast.sam" | tr -d ' ')"

for kind in 2a sa msa; do
	./minibwa bench -b "$kind" -n "$BENCH_ITERATIONS" "$INDEX_FILE" >/dev/null 2>"$TMP/bench-$kind.err"
	append_bwt_timing "bwt_${kind}_${BENCH_ITERATIONS}" "$TMP/bench-$kind.err" "index=$INDEX_FILE"
done

printf "wrote %s\n" "$OUT"
