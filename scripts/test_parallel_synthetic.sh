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

if ! TELSEQ_BENCH_OUT="$results_dir" \
    TELSEQ_STOP_ON_MISMATCH=0 \
    "$repo_dir/scripts/compare_and_benchmark.sh" \
        "$stock_telseq" \
        "$new_telseq" \
        "$bam" \
        1 2 22 44; then
    summary="$results_dir/summary.tsv"
    run_status() {
        awk -F '\t' -v label="$1" '$1 == label { print $2; exit }' \
            "$summary"
    }

    if [ "$(run_status stock)" != "0" ]; then
        echo "FAIL[20]: stock TelSeq execution failed" >&2
        echo "Captured stderr:" >&2
        sed -n '1,200p' "$results_dir/stock.stderr" >&2
        exit 20
    fi

    diagnostic_code=21
    for thread_count in 1 2 22 44; do
        label="threads-$thread_count"
        if [ "$(run_status "$label")" != "0" ]; then
            echo "FAIL[$diagnostic_code]: $label execution failed" >&2
            echo "Captured stderr:" >&2
            sed -n '1,200p' "$results_dir/$label.stderr" >&2
            exit "$diagnostic_code"
        fi
        diagnostic_code=$((diagnostic_code + 1))
    done

    diagnostic_code=31
    for thread_count in 1 2 22 44; do
        label="threads-$thread_count"
        if ! cmp -s \
            "$results_dir/stock.stdout" \
            "$results_dir/$label.stdout"; then
            echo "FAIL[$diagnostic_code]: $label output differs from stock" >&2
            exit "$diagnostic_code"
        fi
        diagnostic_code=$((diagnostic_code + 1))
    done

    echo "FAIL[39]: compatibility benchmark returned an unexplained failure" >&2
    exit 39
fi

if ! grep -q \
    "using 21 mapped-reference workers and 1 HTSlib compatibility scanner across 64 reference tasks" \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: -t 22 did not start the expected worker layout" >&2
    exit 11
fi

if ! grep -q \
    "using 43 mapped-reference workers and 1 HTSlib compatibility scanner across 64 reference tasks" \
    "$results_dir/threads-44.stderr"; then
    echo "FAIL: -t 44 did not start the expected worker layout" >&2
    exit 12
fi

expected_fast_scan="[scan] HTSlib compatibility scanner fetched 8 indexed BAM records: 8 no-coordinate, 0 final-reference fallback; full sequential scans: 0"
if ! grep -Fqx \
    "$expected_fast_scan" \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: compatibility scanner did not use direct indexed tail access" >&2
    echo "Expected log: $expected_fast_scan" >&2
    exit 13
fi

if ! grep -Fqx \
    "[scan] parallel workers processed 1130 reads after filters" \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: indexed tail scan did not retain the legacy final record" >&2
    exit 14
fi

observed_counts=$(awk -F '\t' '$1 == "rg1" {
    print $4 "\t" $5 "\t" $6
}' "$results_dir/stock.stdout")

if [ "$observed_counts" != "$(printf '1130\t1120\t3')" ]; then
    echo "FAIL: expected Total=1130, Mapped=1120, Duplicates=3" >&2
    echo "Observed: $observed_counts" >&2
    exit 15
fi

# Exercise the compatibility fallback used when a coordinate-sorted BAM has
# no no-coordinate tail. The highest populated reference contains 20 records,
# so the scanner must fetch only those records to retain the stock EOF
# contribution; it must not fall back to a whole-file pass.
no_tail_bam="$test_dir/no-tail-fixture.bam"
"$fixture_generator" "$no_tail_bam" 20 100 0

if ! "$stock_telseq" "$no_tail_bam" \
    >"$test_dir/no-tail-stock.stdout" \
    2>"$test_dir/no-tail-stock.stderr"; then
    echo "FAIL[16]: stock TelSeq failed on no-tail fixture" >&2
    exit 16
fi
if ! "$new_telseq" -t 22 "$no_tail_bam" \
    >"$test_dir/no-tail-parallel.stdout" \
    2>"$test_dir/no-tail-parallel.stderr"; then
    echo "FAIL[17]: parallel TelSeq failed on no-tail fixture" >&2
    exit 17
fi

if ! cmp -s \
    "$test_dir/no-tail-stock.stdout" \
    "$test_dir/no-tail-parallel.stdout"; then
    diff -u \
        "$test_dir/no-tail-stock.stdout" \
        "$test_dir/no-tail-parallel.stdout" \
        >"$test_dir/no-tail.diff" || true
    echo "FAIL: no-tail fallback output differs from stock TelSeq" >&2
    echo "Diff: $test_dir/no-tail.diff" >&2
    exit 18
fi

expected_fallback_scan="[scan] HTSlib compatibility scanner fetched 20 indexed BAM records: 0 no-coordinate, 20 final-reference fallback; full sequential scans: 0"
if ! grep -Fqx \
    "$expected_fallback_scan" \
    "$test_dir/no-tail-parallel.stderr"; then
    echo "FAIL: no-tail compatibility path read more than the final reference" >&2
    echo "Expected log: $expected_fallback_scan" >&2
    exit 19
fi

echo
echo "PASS: synthetic parallel test completed"
echo "Expected legacy counts: Total=1130 Mapped=1120 Duplicates=3"
echo "Fast tail access: 8 indexed records fetched; 0 full sequential scans"
echo "No-tail fallback: 20 final-reference records fetched; stock output matched"
echo "Artifacts: $test_dir"
