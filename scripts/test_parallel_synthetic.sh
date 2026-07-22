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

    diagnose_parallel_stderr() {
        stderr_file=$1

        if grep -Fq "HTSlib index does not support direct no-coordinate access" \
            "$stderr_file"; then
            echo 40
        elif grep -Fq "HTSlib could not load a BAI or CSI index" \
            "$stderr_file"; then
            echo 41
        elif grep -Fq "HTSlib failed while reading no-coordinate BAM records" \
            "$stderr_file"; then
            echo 42
        elif grep -Fq "HTSlib could not open BAM" "$stderr_file" ||
             grep -Fq "HTSlib could not read BAM header" "$stderr_file" ||
             grep -Fq "HTSlib could not allocate BAM records" "$stderr_file"; then
            echo 43
        elif grep -Fq "could not copy RG tag from HTSlib alignment" \
            "$stderr_file" ||
             grep -Fq "HTSlib found a non-string RG tag" "$stderr_file"; then
            echo 44
        elif grep -Eq "Error in worker [0-9]+: could not open BAM" \
            "$stderr_file"; then
            echo 45
        elif grep -Eq "Error in worker [0-9]+: could not locate BAM index" \
            "$stderr_file"; then
            echo 46
        elif grep -Eq "Error in worker [0-9]+: could not seek to reference" \
            "$stderr_file"; then
            echo 47
        elif grep -Fq "Command terminated by signal" "$stderr_file"; then
            echo 48
        else
            echo 49
        fi
    }

    if [ "$(run_status stock)" != "0" ]; then
        echo "FAIL[20]: stock TelSeq execution failed" >&2
        echo "Captured stderr:" >&2
        sed -n '1,200p' "$results_dir/stock.stderr" >&2
        exit 20
    fi

    for thread_count in 1 2 22 44; do
        label="threads-$thread_count"
        if [ "$(run_status "$label")" != "0" ]; then
            parallel_diagnostic=$(diagnose_parallel_stderr \
                "$results_dir/$label.stderr")
            echo "FAIL[$parallel_diagnostic]: $label execution failed" >&2
            echo "Captured stderr:" >&2
            sed -n '1,200p' "$results_dir/$label.stderr" >&2
            exit "$parallel_diagnostic"
        fi
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
    "using 21 mapped-reference workers and 1 HTSlib compatibility scanner across 64 mapped-reference tasks from 64 references; window size 25000000 bp" \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: -t 22 did not start the expected worker layout" >&2
    exit 11
fi

if ! grep -q \
    "using 43 mapped-reference workers and 1 HTSlib compatibility scanner across 64 mapped-reference tasks from 64 references; window size 25000000 bp" \
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

if grep -q '^\[reference-profile\]' \
    "$results_dir/threads-22.stderr"; then
    echo "FAIL: reference profiling was enabled by default" >&2
    exit 24
fi

profile_stdout="$test_dir/profile.stdout"
profile_stderr="$test_dir/profile.stderr"
if ! "$new_telseq" --profile-references -t 22 "$bam" \
    >"$profile_stdout" \
    2>"$profile_stderr"; then
    echo "FAIL[25]: profiled parallel execution failed" >&2
    sed -n '1,200p' "$profile_stderr" >&2
    exit 25
fi

if ! cmp -s "$results_dir/threads-22.stdout" "$profile_stdout"; then
    echo "FAIL: reference profiling changed standard output" >&2
    exit 26
fi

if ! awk -F '\t' '
    $1 == "[reference-profile]" && $2 == "task" {
        headers += 1
        next
    }
    $1 == "[reference-profile]" {
        rows += 1
        if (NF != 13 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ ||
            $3 < 1 || $3 > 21 || $4 !~ /^[0-9]+$/ ||
            $6 !~ /^[0-9]+$/ || $7 !~ /^[0-9]+$/ ||
            $8 !~ /^[0-9]+$/ || $9 !~ /^[0-9]+$/ ||
            $10 !~ /^[0-9]+$/ || $11 < 0 || $12 < $11 || $13 < 0) {
            invalid = 1
        }
        seen_reference[$4] += 1
    }
    END {
        if (headers != 1 || rows != 64 || invalid) {
            exit 1
        }
        for (reference in seen_reference) {
            if (seen_reference[reference] != 1) {
                exit 1
            }
        }
    }
' "$profile_stderr"; then
    echo "FAIL: per-reference profile is incomplete or malformed" >&2
    sed -n '/^\[reference-profile\]/p' "$profile_stderr" >&2
    exit 27
fi

# Split one long synthetic reference into three windows. Records ending at,
# spanning, starting at, and starting immediately after a boundary prove that
# overlap queries are reduced to exact start-coordinate ownership.
window_bam="$test_dir/window-boundary-fixture.bam"
window_stock_stdout="$test_dir/window-boundary-stock.stdout"
window_parallel_stdout="$test_dir/window-boundary-parallel.stdout"
window_parallel_stderr="$test_dir/window-boundary-parallel.stderr"
"$fixture_generator" "$window_bam" 20 100 8 1

if ! "$stock_telseq" "$window_bam" \
    >"$window_stock_stdout" \
    2>"$test_dir/window-boundary-stock.stderr"; then
    echo "FAIL[28]: stock TelSeq failed on window-boundary fixture" >&2
    exit 28
fi
if ! "$new_telseq" --profile-references -t 22 "$window_bam" \
    >"$window_parallel_stdout" \
    2>"$window_parallel_stderr"; then
    echo "FAIL[29]: windowed parallel execution failed" >&2
    sed -n '1,200p' "$window_parallel_stderr" >&2
    exit 29
fi

if ! cmp -s "$window_stock_stdout" "$window_parallel_stdout"; then
    diff -u "$window_stock_stdout" "$window_parallel_stdout" \
        >"$test_dir/window-boundary.diff" || true
    echo "FAIL: window-boundary output differs from stock TelSeq" >&2
    echo "Diff: $test_dir/window-boundary.diff" >&2
    exit 30
fi

if ! awk -F '\t' '
    $1 == "[reference-profile]" && $2 ~ /^[0-9]+$/ {
        rows += 1
        if ($5 == "contig0") {
            contig_tasks += 1
            if ($7 == 0 && $8 == 25000000 && $9 == 23) first = 1
            if ($7 == 25000000 && $8 == 50000000 && $9 == 3) second = 1
            if ($7 == 50000000 && $8 == 50001100 && $9 == 1) third = 1
        }
    }
    END {
        if (rows != 66 || contig_tasks != 3 ||
            !first || !second || !third) {
            exit 1
        }
    }
' "$window_parallel_stderr"; then
    echo "FAIL: window ownership duplicated or omitted boundary records" >&2
    sed -n '/^\[reference-profile\]/p' "$window_parallel_stderr" >&2
    exit 31
fi

window_counts=$(awk -F '\t' '$1 == "rg1" {
    print $4 "\t" $5 "\t" $6
}' "$window_stock_stdout")
if [ "$window_counts" != "$(printf '1136\t1125\t3')" ]; then
    echo "FAIL: unexpected counts in window-boundary fixture" >&2
    echo "Observed: $window_counts" >&2
    exit 32
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
echo "Reference profile: 64 timed tasks; standard output unchanged"
echo "Window ownership: boundary-spanning records counted exactly once"
echo "Artifacts: $test_dir"
