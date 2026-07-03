#!/usr/bin/env sh
set -eu

mkdir -p build
: "${CFLAGS:=-std=c99 -g -Wall -O0}"
cc $CFLAGS -Iinclude -Isrc tests/test-par.c libminibwa.a -lz -lm -lpthread -o build/test-par
./build/test-par
