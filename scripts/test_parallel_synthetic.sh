#!/usr/bin/env bash

set -eu

usage() {
    echo "Usage: $0 STOCK_TELSEQ NEW_TELSEQ [FIXTURE_GENERATOR]" >&2
}

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    usage
    exit 2
fi

stock_telseq=$1
new_telseq=$2

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "$script_dir/.." && pwd)
fixture_generator=${3:-"$repo_dir/src/Test/generate_parallel_fixture"}

for executable in "$stock_telseq" "$new_telseq" "$fixture_generator"; do
    if [ ! -x "$executable" ]; then
        echo "Error: executable not found: $executable" >&2
        exit 2
    fi
done

test_dir=$(mktemp -d "${TMPDIR:-/tmp}/telseq-parallel-test.XXXXXX")
bam="$test_dir/parallel-fixture.bam"
results_dir="$test_dir/results"

"$fixture_generator" "$bam"

TELSEQ_BENCH_OUT="$results_dir" \
    "$repo_dir/scripts/compare_and_benchmark.sh" \
    "$stock_telseq" \
    "$new_telseq" \
    "$bam" \
    1 2 22 44

if ! grep -q \
    "using 21 mapped-reference workers and 1 no-coordinate scanner across 64 reference tasks" \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: -t 22 did not start the expected worker layout" >&2
    exit 1
fi

if ! grep -q \
    "using 43 mapped-reference workers and 1 no-coordinate scanner across 64 reference tasks" \
    "$results_dir/threads-44.stderr"; then
    echo "FAIL: -t 44 did not start the expected worker layout" >&2
    exit 1
fi

observed_counts=$(awk -F '\t' '$1 == "rg1" {
    print $4 "\t" $5 "\t" $6
}' "$results_dir/stock.stdout")

if [ "$observed_counts" != "$(printf '1130\t1120\t3')" ]; then
    echo "FAIL: expected Total=1130, Mapped=1120, Duplicates=3" >&2
    echo "Observed: $observed_counts" >&2
    exit 1
fi

echo
echo "PASS: synthetic parallel test completed"
echo "Expected legacy counts: Total=1130 Mapped=1120 Duplicates=3"
echo "Artifacts: $test_dir"
