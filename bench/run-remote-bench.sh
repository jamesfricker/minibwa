#!/usr/bin/env bash
set -euo pipefail

usage() {
	cat <<'EOF'
usage: bench/run-remote-bench.sh <ssh-host>

Copies this source tree to an SSH host, builds it, and optionally runs the
minibwa bench subcommand against a caller-provided remote index.

Environment:
  REMOTE_BASE       Remote work parent directory (default: ${HOME}/minibwa-bench)
  REMOTE_WORK       Full remote work directory (default: $REMOTE_BASE/run-<timestamp>)
  BENCH_INDEX       Remote .mbw index to benchmark. If unset, only build is run.
  BENCH_ITERATIONS  Iterations for each bench mode (default: 1000000)
  LOCAL_LOG         Local log path (default: .context/remote-bench-<timestamp>.log)
EOF
}

ROOT=$(cd "$(dirname "$0")/.." && pwd)
REMOTE=${1:-${REMOTE:-}}
if [ -z "$REMOTE" ]; then
	usage >&2
	exit 2
fi

STAMP=$(date -u +"%Y%m%dT%H%M%SZ")
REMOTE_BASE=${REMOTE_BASE:-'${HOME}/minibwa-bench'}
REMOTE_WORK=${REMOTE_WORK:-"$REMOTE_BASE/run-$STAMP"}
BENCH_INDEX=${BENCH_INDEX:-}
BENCH_ITERATIONS=${BENCH_ITERATIONS:-1000000}
LOCAL_LOG=${LOCAL_LOG:-"$ROOT/.context/remote-bench-$STAMP.log"}

cd "$ROOT"
mkdir -p "$ROOT/.context"

COPYFILE_DISABLE=1 tar --no-xattrs --exclude .git --exclude .context --exclude '*.o' --exclude '*.a' --exclude minibwa -czf - . \
	| ssh "$REMOTE" "set -euo pipefail; mkdir -p \"$REMOTE_WORK/src\"; tar -xzf - -C \"$REMOTE_WORK/src\""

ssh "$REMOTE" "set -euo pipefail
	cd \"$REMOTE_WORK/src\"
	make clean >/dev/null
	make >/dev/null
	if [ -n \"$BENCH_INDEX\" ]; then
		if [ ! -f \"$BENCH_INDEX\" ]; then
			printf 'benchmark index not found: %s\n' \"$BENCH_INDEX\" >&2
			exit 1
		fi
		for kind in 2a sa msa; do
			/usr/bin/time -f \"bench_\${kind}_elapsed_seconds\t%e\" \
				\"$REMOTE_WORK/src/minibwa\" bench -b \"\$kind\" -n \"$BENCH_ITERATIONS\" \"$BENCH_INDEX\" >/dev/null
		done
	fi
" 2>&1 | tee "$LOCAL_LOG"

printf "remote workdir: %s:%s\nlocal log: %s\n" "$REMOTE" "$REMOTE_WORK" "$LOCAL_LOG"
