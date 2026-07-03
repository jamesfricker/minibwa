#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${MINIBWA:-$root/minibwa}"
tmp="${TMPDIR:-/tmp}/minibwa-human-ambiguous.$$"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp"

cat > "$tmp/ref.fa" <<'FA'
>ambiguous
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
FA

"$bin" index "$tmp/ref.fa" "$tmp/default" >/dev/null 2>"$tmp/default.index.log"
"$bin" index --human "$tmp/ref.fa" "$tmp/human" >/dev/null 2>"$tmp/human.index.log"

python3 - "$tmp/default.l2b" "$tmp/read.fa" <<'PY'
import struct
import sys

l2b_path, read_path = sys.argv[1:3]
with open(l2b_path, "rb") as fh:
    magic = fh.read(4)
    if magic != b"L2B\1":
        raise SystemExit("bad l2b magic")
    fh.read(4)
    n_ctg, tot_len, n_ambi, n_mask, len_name, len_comm, n_pac = struct.unpack("<7Q", fh.read(56))
    fh.read(8 * n_ctg)
    fh.read(16 * n_ambi)
    fh.read(16 * n_mask)
    pac = fh.read(8 * n_pac)

bases = "ACGT"
seq = []
for i in range(tot_len):
    word = struct.unpack_from("<Q", pac, (i >> 5) * 8)[0]
    seq.append(bases[(word >> ((i & 31) * 2)) & 3])

with open(read_path, "w") as out:
    out.write(">read_from_randomized_ns\n")
    out.write("".join(seq) + "\n")
PY

"$bin" map -f -u --chain-only "$tmp/default" "$tmp/read.fa" > "$tmp/default.paf"
"$bin" map -f -u --chain-only "$tmp/human" "$tmp/read.fa" > "$tmp/human.paf"

if [[ ! -s "$tmp/default.paf" ]]; then
	echo "expected default index to map through randomized N bases" >&2
	exit 1
fi

if [[ -s "$tmp/human.paf" ]]; then
	echo "expected human index not to map through randomized N bases" >&2
	cat "$tmp/human.paf" >&2
	exit 1
fi
