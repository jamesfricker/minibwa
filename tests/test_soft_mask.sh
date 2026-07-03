#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

make -s minibwa libminibwa.a

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/ref.fa" <<'FA'
>chr1
acacacacacacacacacacacacacacacacacacacacacacacacacacacacacacacacGTGTTGCATGTCAGTACCGTAGCTAGGCTTACGATCGTACGATGCTAGCATCGATCGTACccgtaacgttgacctgatcgtacgatgctagtcagtacgatcgtagctagtcagtcgatgacctgatcgtagctagttcgatcgatgctagtcgatcgtac
FA

cat > "$tmp/reads.fa" <<'FA'
>upper_unique
GTGTTGCATGTCAGTACCGTAGCTAGGCTTACGATCGTACGATGCTAGCATCG
>terminal_masked_unique
CCGTAACGTTGACCTGATCGTACGATGCTAGTCAGTACGATCGTAGCTAGTCAG
>masked_repeat
ACACACACACACACACACACACACACACACACACACACAC
FA

./minibwa index "$tmp/ref.fa" "$tmp/ref" >/dev/null 2>&1

cat > "$tmp/probe.c" <<'C'
#include <stdio.h>
#include "l2bit.h"

int main(int argc, char **argv)
{
	l2b_t *l2b;
	int failed = 0;
	if (argc != 2) return 2;
	l2b = l2b_load(argv[1]);
	if (l2b == 0) return 2;
	if (l2b_mask_overlap(l2b, 0, 0, 40) != 40) failed = 1;
	if (l2b_mask_overlap(l2b, 0, 64, 117) != 0) failed = 1;
	if (l2b_mask_overlap(l2b, 0, 124, 178) != 54) failed = 1;
	l2b_destroy(l2b);
	return failed;
}
C

${CC:-cc} -std=c99 -Iinclude -Isrc "$tmp/probe.c" libminibwa.a -lpthread -lz -lm -o "$tmp/probe"
"$tmp/probe" "$tmp/ref.l2b"

./minibwa map -f -P -x sr -t 1 --dbg-anchor "$tmp/ref" "$tmp/reads.fa" \
	> "$tmp/out.paf" 2> "$tmp/debug.log"

upper_mapq=$(awk '$1 == "upper_unique" { print $12 }' "$tmp/out.paf")
terminal_mapq=$(awk '$1 == "terminal_masked_unique" { print $12 }' "$tmp/out.paf")
repeat_target=$(awk '$1 == "masked_repeat" { print $6 }' "$tmp/out.paf")

test "$upper_mapq" = 60
test "$terminal_mapq" = 10
test "$repeat_target" = "*"
if grep -q $'^AC\tmasked_repeat\t' "$tmp/debug.log"; then
	echo "masked_repeat unexpectedly produced anchors" >&2
	exit 1
fi

echo "soft-mask test passed"
