#!/bin/sh
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

make -s minibwa

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat >"$tmp/ref.fa" <<'FA'
>chr1
ACGTTGCAAGTCGATCGTACCGATGCTAGTACGATCGATGGTACCGTACGTTAGCCGATACGTTGACCGTAGCTTACGATGGCATCGTACCGTTAACGATCGTAGCTAGGCTAACCGTATGCATC
FA

cat >"$tmp/reads.fa" <<'FA'
>read_alpha
ACGTTGCAAGTCGATCGTACCGATGCTAGTACGAT
>read_beta
TTGACCGTAGCTTACGATGGCATCGTACCGTTAA
FA

./minibwa index "$tmp/ref.fa" "$tmp/ref" >/dev/null 2>"$tmp/index.log"
./minibwa map -t1 -x sr "$tmp/ref" "$tmp/reads.fa" >"$tmp/default.sam" 2>"$tmp/default.log"
./minibwa map -t1 -x sr -P "$tmp/ref" "$tmp/reads.fa" >"$tmp/P.sam" 2>"$tmp/P.log"
./minibwa map -t1 -x sr --single-end "$tmp/ref" "$tmp/reads.fa" >"$tmp/single.sam" 2>"$tmp/single.log"

./minibwa map --help >"$tmp/help.txt" 2>"$tmp/help.log"
grep -q -- "--single-end" "$tmp/help.txt"

grep -v '^@PG' "$tmp/default.sam" >"$tmp/default.body"
grep -v '^@PG' "$tmp/P.sam" >"$tmp/P.body"
grep -v '^@PG' "$tmp/single.sam" >"$tmp/single.body"

cmp "$tmp/default.body" "$tmp/P.body"
cmp "$tmp/P.body" "$tmp/single.body"

echo "single-end fast-path option test passed"
