#!/usr/bin/env bash

set -eu

usage() {
    echo "Usage: $0 TELSEQ FIXTURE_GENERATOR" >&2
}

if [ "$#" -ne 2 ]; then
    usage
    exit 2
fi

telseq=$1
fixture_generator=$2
reads_per_reference=${TELSEQ_SCALING_READS_PER_REFERENCE:-4000}
read_length=${TELSEQ_SCALING_READ_LENGTH:-1000}

for executable in "$telseq" "$fixture_generator"; do
    if [ ! -x "$executable" ]; then
        echo "Error: executable not found: $executable" >&2
        exit 2
    fi
done

if ! command -v nproc >/dev/null 2>&1 || [ "$(nproc)" -lt 2 ]; then
    echo "SKIP: scaling regression requires at least two logical CPUs" >&2
    exit 0
fi

test_dir=$(mktemp -d "${TMPDIR:-/tmp}/telseq-scaling-test.XXXXXX")
bam="$test_dir/scaling-fixture.bam"

"$fixture_generator" "$bam" "$reads_per_reference" "$read_length"

run_timed() {
    threads=$1
    run=$2
    stdout="$test_dir/threads-$threads-run-$run.stdout"
    stderr="$test_dir/threads-$threads-run-$run.stderr"

    start_ns=$(date +%s%N)
    if ! "$telseq" -t "$threads" -u "$bam" >"$stdout" 2>"$stderr"; then
        echo "FAIL[50]: -t $threads execution failed" >&2
        sed -n '1,200p' "$stderr" >&2
        return 50
    fi
    end_ns=$(date +%s%N)

    awk -v start="$start_ns" -v end="$end_ns" \
        'BEGIN { printf "%.6f\n", (end - start) / 1000000000 }'
}

# Warm both layouts before measuring so page-cache state is not assigned to
# only one thread count.
run_timed 2 warm >/dev/null
run_timed 4 warm >/dev/null

: >"$test_dir/threads-2.times"
: >"$test_dir/threads-4.times"

for run in 1 2 3; do
    run_timed 2 "$run" >>"$test_dir/threads-2.times"
    run_timed 4 "$run" >>"$test_dir/threads-4.times"
done

if ! cmp -s \
    "$test_dir/threads-2-run-1.stdout" \
    "$test_dir/threads-4-run-1.stdout"; then
    echo "FAIL[51]: -t 2 and -t 4 outputs differ" >&2
    exit 51
fi

median_seconds() {
    sort -n "$1" | awk 'NR == 2 { print; exit }'
}

threads_2_median=$(median_seconds "$test_dir/threads-2.times")
threads_4_median=$(median_seconds "$test_dir/threads-4.times")

if ! awk -v slower="$threads_2_median" -v faster="$threads_4_median" \
    'BEGIN { exit !(faster < slower) }'; then
    echo "FAIL: additional mapped-reference workers did not reduce median wall time" >&2
    echo "-t 2 median: $threads_2_median seconds" >&2
    echo "-t 4 median: $threads_4_median seconds" >&2
    echo "Artifacts: $test_dir" >&2
    if awk -v slower="$threads_2_median" -v faster="$threads_4_median" \
        'BEGIN { exit !(faster <= slower * 1.05) }'; then
        exit 52
    elif awk -v slower="$threads_2_median" -v faster="$threads_4_median" \
        'BEGIN { exit !(faster <= slower * 1.20) }'; then
        exit 53
    else
        exit 54
    fi
fi

speedup=$(awk -v slower="$threads_2_median" -v faster="$threads_4_median" \
    'BEGIN { printf "%.2f", slower / faster }')

echo "PASS: additional threads reduced synthetic workload wall time"
echo "-t 2 median: $threads_2_median seconds"
echo "-t 4 median: $threads_4_median seconds"
echo "Observed speedup: ${speedup}x"
echo "Artifacts: $test_dir"
