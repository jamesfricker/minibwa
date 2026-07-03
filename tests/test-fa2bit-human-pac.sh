#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bin="${MINIBWA:-$root/minibwa}"
tmp="${TMPDIR:-/tmp}/minibwa-fa2bit-human-pac.$$"
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp"

cat > "$tmp/ref.fa" <<'FA'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGT
FA

# 1. fa2bit --human combined with -p must be rejected: exit 1, error on stderr,
#    and no output file written (BWA pac format cannot store the human flag).
set +e
"$bin" fa2bit --human -p "$tmp/ref.fa" "$tmp/out.pac" >"$tmp/rej.out" 2>"$tmp/rej.err"
rc=$?
set -e
if [[ "$rc" -ne 1 ]]; then
	echo "expected fa2bit --human -p to exit 1, got $rc" >&2
	exit 1
fi
if ! grep -q "not supported with -p" "$tmp/rej.err"; then
	echo "expected fa2bit --human -p to print a rejection message on stderr" >&2
	cat "$tmp/rej.err" >&2
	exit 1
fi
if [[ -e "$tmp/out.pac" ]]; then
	echo "expected no output file when fa2bit --human -p is rejected" >&2
	exit 1
fi

# 2. fa2bit --human on its own (l2b output) must still succeed.
"$bin" fa2bit --human "$tmp/ref.fa" "$tmp/out.l2b" >/dev/null 2>&1
if [[ ! -s "$tmp/out.l2b" ]]; then
	echo "expected fa2bit --human to write an l2b file" >&2
	exit 1
fi

# 3. fa2bit -p on its own (pac output) must still succeed.
"$bin" fa2bit -p "$tmp/ref.fa" "$tmp/ok.pac" >/dev/null 2>&1
if [[ ! -s "$tmp/ok.pac" ]]; then
	echo "expected fa2bit -p to write a pac file" >&2
	exit 1
fi

# 4. Exit codes from subcommands must propagate through top-level main().
#    Before the fix, main() returned 0 unconditionally and masked failures.
set +e
"$bin" fa2bit --human -p "$tmp/ref.fa" "$tmp/out2.pac" >/dev/null 2>&1
top_rc=$?
set -e
if [[ "$top_rc" -ne 1 ]]; then
	echo "expected top-level exit code 1 to propagate from fa2bit, got $top_rc" >&2
	exit 1
fi

echo "fa2bit --human/-p rejection and exit-code propagation test passed"
