#!/usr/bin/env bash

set -u

usage() {
    echo "Usage: $0 STOCK_TELSEQ NEW_TELSEQ SAMPLE_BAM [THREADS ...]" >&2
    echo "Example: $0 ./stock-telseq ./src/Telseq/telseq sample.bam 1 2 4 8" >&2
}

if [ "$#" -lt 3 ]; then
    usage
    exit 2
fi

stock_telseq=$1
new_telseq=$2
bam=$3
shift 3

if [ "$#" -eq 0 ]; then
    set -- 1 2 4 8 16 24
fi

for required_file in "$stock_telseq" "$new_telseq" "$bam"; do
    if [ ! -e "$required_file" ]; then
        echo "Error: not found: $required_file" >&2
        exit 2
    fi
done

if [ ! -x "$stock_telseq" ]; then
    echo "Error: stock TelSeq is not executable: $stock_telseq" >&2
    exit 2
fi
if [ ! -x "$new_telseq" ]; then
    echo "Error: new TelSeq is not executable: $new_telseq" >&2
    exit 2
fi

case "$bam" in
    *.bam)
        bam_without_suffix=${bam%.bam}
        ;;
    *)
        bam_without_suffix=$bam
        ;;
esac

if [ ! -e "${bam}.bai" ] && [ ! -e "${bam_without_suffix}.bai" ] &&
   [ ! -e "${bam}.bti" ]; then
    echo "Warning: no .bai or .bti was found next to $bam." >&2
    echo "Parallel runs are expected to fail until an index is available." >&2
fi

timestamp=$(date +%Y%m%d-%H%M%S)
output_dir=${TELSEQ_BENCH_OUT:-"telseq-benchmark-${timestamp}"}
mkdir -p "$output_dir"

summary="$output_dir/summary.tsv"
environment_log="$output_dir/environment.txt"

{
    echo "date: $(date)"
    echo "host: $(hostname 2>/dev/null || echo unknown)"
    echo "system: $(uname -a)"
    echo "stock: $stock_telseq"
    echo "new: $new_telseq"
    echo "bam: $bam"
    echo "threads: $*"
    echo
    echo "new TelSeq version:"
    "$new_telseq" --version 2>&1 || true
    if command -v samtools >/dev/null 2>&1; then
        echo
        echo "samtools version:"
        samtools --version 2>&1 | sed -n '1,2p'
        echo
        echo "samtools quickcheck:"
        samtools quickcheck -v "$bam" 2>&1
        echo "samtools physical alignment count:"
        samtools view -c "$bam" 2>&1
    fi
} > "$environment_log"

checksum_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

run_timed() {
    label=$1
    shift

    stdout_file="$output_dir/${label}.stdout"
    stderr_file="$output_dir/${label}.stderr"

    /usr/bin/time -p "$@" > "$stdout_file" 2> "$stderr_file"
    run_status=$?

    real_time=$(awk '$1 == "real" {value=$2} END {print value}' "$stderr_file")
    user_time=$(awk '$1 == "user" {value=$2} END {print value}' "$stderr_file")
    sys_time=$(awk '$1 == "sys" {value=$2} END {print value}' "$stderr_file")
    checksum=$(checksum_file "$stdout_file")

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$run_status" "${real_time:-NA}" "${user_time:-NA}" \
        "${sys_time:-NA}" "$checksum" >> "$summary"

    return "$run_status"
}

printf 'run\texit_status\treal_seconds\tuser_seconds\tsys_seconds\tstdout_sha256\n' \
    > "$summary"

echo "Running stock TelSeq..."
if ! run_timed stock "$stock_telseq" "$bam"; then
    echo "Error: stock TelSeq failed. See $output_dir/stock.stderr" >&2
    echo "Results: $output_dir"
    exit 1
fi

stop_on_mismatch=${TELSEQ_STOP_ON_MISMATCH:-1}
had_failure=0

for threads in "$@"; do
    case "$threads" in
        ''|*[!0-9]*)
            echo "Error: invalid thread count: $threads" >&2
            exit 2
            ;;
    esac

    label="threads-${threads}"
    echo "Running new TelSeq with -t $threads..."
    if ! run_timed "$label" "$new_telseq" -t "$threads" "$bam"; then
        echo "FAIL: -t $threads exited unsuccessfully; see ${label}.stderr" >&2
        had_failure=1
        if [ "$stop_on_mismatch" != "0" ]; then
            break
        fi
        continue
    fi

    if cmp -s "$output_dir/stock.stdout" "$output_dir/${label}.stdout"; then
        echo "PASS: -t $threads output is byte-identical to stock"
    else
        echo "FAIL: -t $threads output differs from stock" >&2
        diff -u "$output_dir/stock.stdout" "$output_dir/${label}.stdout" \
            > "$output_dir/${label}.diff" || true
        echo "Diff: $output_dir/${label}.diff" >&2
        had_failure=1
        if [ "$stop_on_mismatch" != "0" ]; then
            break
        fi
    fi
done

echo
echo "Summary:"
column -t -s $'\t' "$summary" 2>/dev/null || cat "$summary"
echo
echo "Results: $output_dir"

exit "$had_failure"
