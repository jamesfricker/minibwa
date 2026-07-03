#!/bin/sh
set -eu

MINIBWA=${MINIBWA:-./minibwa}
if [ ! -x "$MINIBWA" ]; then
	echo "missing executable: $MINIBWA" >&2
	exit 1
fi

tmp=$(mktemp -d "${TMPDIR:-/tmp}/minibwa-human-tags.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

cat > "$tmp/ref.fa" <<'FA'
>chr1
ACGTTGCAAGTCGATCGTACCGATGCTAGTACGATCGATGGTACCGTACGTTAGCCGATACGT
>chr6_GL000250v2_alt_HLA
TTGACCGTAGCTTACGATGGCATCGTACCGTTAACGATCGTAGCTAGGCTAACCGTATGCATC
>chr2_KI270715v1_random
GATCCGTAACGTTAGCATCGGATACGTCAGTTAACCGGTTACGATCGTAGCTTACCGATGCTA
>chrUn_GL000220v1
CCGTAGTACGATTCGATGCAACCGTAGGCTAACGTTAGCATCGATACGGCATTCGATCGTACC
>hs37d5
TATGCGATACCGTTAAGCTAGCATCGTACGATGGCCTAACGATCGTATTCGGCATAGCTACGA
>chr3
GGCATTCGATCGTACCGTAAttgcaacgtagctacgatggcatcgtaccgttAACCGGTTAGCATCGATGCA
FA

cat > "$tmp/reads.fa" <<'FA'
>primary
ACGTTGCAAGTCGATCGTACCGATGCTAGTACGAT
>hla
TTGACCGTAGCTTACGATGGCATCGTACCGTTAA
>random
GATCCGTAACGTTAGCATCGGATACGTCAGTTAA
>unplaced
CCGTAGTACGATTCGATGCAACCGTAGGCTAACG
>decoy
TATGCGATACCGTTAAGCTAGCATCGTACGATGG
>low
TTGCAACGTAGCTACGATGGCATCGTACCGTT
FA

"$MINIBWA" index "$tmp/ref.fa" "$tmp/ref" >/dev/null 2>"$tmp/index.log"
"$MINIBWA" map --adap=no --pe=no --human-tags "$tmp/ref" "$tmp/reads.fa" > "$tmp/out.sam" 2>"$tmp/sam.log"
"$MINIBWA" map --adap=no --pe=no -f --human-tags "$tmp/ref" "$tmp/reads.fa" > "$tmp/out.paf" 2>"$tmp/paf.log"
"$MINIBWA" map --adap=no --pe=no "$tmp/ref" "$tmp/reads.fa" > "$tmp/default.sam" 2>"$tmp/default.log"

assert_record_has() {
	file=$1
	name=$2
	text=$3
	if ! awk -v name="$name" -v text="$text" '$1 == name && index($0, text) { found = 1 } END { exit found ? 0 : 1 }' "$file"; then
		echo "expected $name record in $file to contain $text" >&2
		echo "--- $file ---" >&2
		cat "$file" >&2
		exit 1
	fi
}

assert_record_has "$tmp/out.sam" primary "zc:Z:primary"
assert_record_has "$tmp/out.sam" primary "zh:Z:high_confidence"
assert_record_has "$tmp/out.sam" hla "zc:Z:hla"
assert_record_has "$tmp/out.sam" random "zc:Z:random"
assert_record_has "$tmp/out.sam" unplaced "zc:Z:unplaced"
assert_record_has "$tmp/out.sam" decoy "zc:Z:decoy"
assert_record_has "$tmp/out.sam" low "zc:Z:primary"
assert_record_has "$tmp/out.sam" low "zm:f:0.000"
assert_record_has "$tmp/out.sam" low "zh:Z:low_mappability"

assert_record_has "$tmp/out.paf" primary "zc:Z:primary"
assert_record_has "$tmp/out.paf" hla "zc:Z:hla"
assert_record_has "$tmp/out.paf" random "zc:Z:random"
assert_record_has "$tmp/out.paf" unplaced "zc:Z:unplaced"
assert_record_has "$tmp/out.paf" decoy "zc:Z:decoy"
assert_record_has "$tmp/out.paf" low "zm:f:0.000"
assert_record_has "$tmp/out.paf" low "zh:Z:low_mappability"

if grep -Eq 'zc:Z:|zm:f:|zh:Z:' "$tmp/default.sam"; then
	echo "human tags were emitted without --human-tags" >&2
	cat "$tmp/default.sam" >&2
	exit 1
fi

if ! grep -q '^@CO	minibwa human-tags:' "$tmp/out.sam"; then
	echo "SAM header is missing human-tags comments" >&2
	cat "$tmp/out.sam" >&2
	exit 1
fi
