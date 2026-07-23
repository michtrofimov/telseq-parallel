#!/usr/bin/env bash

set -u

usage() {
    echo "Usage:" >&2
    echo "  $0 --reference-output STOCK_OUTPUT IMAGE SAMPLE_BAM [THREADS ...] [-- TELSEQ_ARGS ...]" >&2
}

if [ "${1:-}" = "--help" ]; then
    usage
    exit 0
fi

if [ "$#" -lt 4 ] || [ "$1" != "--reference-output" ]; then
    usage
    exit 2
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
output_comparator="$script_dir/compare_telseq_output.sh"
if [ ! -x "$output_comparator" ]; then
    echo "Error: output comparator is not executable: $output_comparator" >&2
    exit 2
fi

reference_output=$2
image=$3
bam=$4
shift 4

thread_grid=()
telseq_args=()

while [ "$#" -gt 0 ]; do
    if [ "$1" = "--" ]; then
        shift
        while [ "$#" -gt 0 ]; do
            telseq_args+=("$1")
            shift
        done
        break
    fi
    thread_grid+=("$1")
    shift
done

if [ "${#thread_grid[@]}" -eq 0 ]; then
    thread_grid=(1 2 4 8 16 24)
fi

for argument in "${telseq_args[@]}"; do
    case "$argument" in
        -t|-t*|--threads|--threads=*)
            echo "Error: pass thread counts before --; do not include $argument in TELSEQ_ARGS" >&2
            exit 2
            ;;
        -f|-f*|--bamlist|--bamlist=*)
            echo "Error: the benchmark controls the input BAM; do not pass $argument" >&2
            exit 2
            ;;
        -o|-o*|--output-dir|--output-dir=*)
            echo "Error: the benchmark captures stdout; do not pass $argument" >&2
            exit 2
            ;;
    esac
done

for thread_count in "${thread_grid[@]}"; do
    case "$thread_count" in
        ''|*[!0-9]*)
            echo "Error: invalid thread count: $thread_count" >&2
            exit 2
            ;;
        0)
            echo "Error: thread count must be at least 1" >&2
            exit 2
            ;;
    esac
done

for required_file in "$reference_output" "$bam"; do
    if [ ! -f "$required_file" ]; then
        echo "Error: not found: $required_file" >&2
        exit 2
    fi
done

docker_bin=${TELSEQ_DOCKER_BIN:-docker}
if ! command -v "$docker_bin" >/dev/null 2>&1; then
    echo "Error: container runtime not found: $docker_bin" >&2
    exit 2
fi
if [ ! -x /usr/bin/time ]; then
    echo "Error: /usr/bin/time is required for benchmarking" >&2
    exit 2
fi
if ! "$docker_bin" info >/dev/null 2>&1; then
    echo "Error: cannot connect to the Docker daemon" >&2
    exit 2
fi

bam_dir=$(CDPATH= cd -- "$(dirname -- "$bam")" && pwd -P)
bam_name=$(basename -- "$bam")
bam="$bam_dir/$bam_name"
# Preserve the host's absolute BAM path inside the container. The inherited
# TelSeq implementation writes the input path to stdout, so mounting it at a
# synthetic path such as /data would create a false output mismatch.
container_bam=$bam

case "$bam" in
    *.bam)
        bam_without_suffix=${bam%.bam}
        ;;
    *)
        bam_without_suffix=$bam
        ;;
esac

if [ ! -e "${bam}.bai" ] && [ ! -e "${bam_without_suffix}.bai" ]; then
    echo "Warning: no BAI was found next to $bam." >&2
    echo "Parallel container runs are expected to fail until an index is available." >&2
fi

if [ "${TELSEQ_DOCKER_PLATFORM+x}" = x ]; then
    docker_platform=$TELSEQ_DOCKER_PLATFORM
else
    docker_platform=linux/amd64
fi

docker_run=("$docker_bin" run --rm)
if [ -n "$docker_platform" ]; then
    docker_run+=(--platform "$docker_platform")
fi
docker_run+=(
    --mount "type=bind,source=$bam_dir,target=$bam_dir,readonly"
    "$image"
)

if ! container_version=$("${docker_run[@]}" --version 2>&1); then
    echo "Error: could not start $image" >&2
    echo "$container_version" >&2
    exit 2
fi

timestamp=$(date +%Y%m%d-%H%M%S)
output_dir=${TELSEQ_BENCH_OUT:-"telseq-docker-benchmark-${timestamp}"}
mkdir -p "$output_dir"

summary="$output_dir/summary.tsv"
environment_log="$output_dir/environment.txt"

{
    echo "date: $(date)"
    echo "host: $(hostname 2>/dev/null || echo unknown)"
    echo "system: $(uname -a)"
    echo "reference mode: file"
    echo "reference: $reference_output"
    echo "container runtime: $docker_bin"
    echo "image: $image"
    echo "platform: ${docker_platform:-runtime default}"
    echo "host BAM: $bam"
    echo "container BAM: $container_bam"
    echo "threads: ${thread_grid[*]}"
    printf 'telseq args:'
    printf ' %q' "${telseq_args[@]}"
    printf '\n'
    echo
    echo "container TelSeq version:"
    echo "$container_version"
    echo
    echo "container image identity:"
    "$docker_bin" image inspect \
        --format 'ID={{.Id}} RepoDigests={{json .RepoDigests}}' \
        "$image" 2>&1 || true
    if command -v samtools >/dev/null 2>&1; then
        echo
        echo "samtools version:"
        samtools --version 2>&1 | sed -n '1,2p'
        echo
        echo "samtools quickcheck:"
        samtools quickcheck -v "$bam" 2>&1
        if [ "${TELSEQ_RUN_BAM_COUNT:-0}" = "1" ]; then
            echo "samtools physical alignment count:"
            samtools view -c "$bam" 2>&1
        else
            echo "samtools physical alignment count: skipped"
            echo "Set TELSEQ_RUN_BAM_COUNT=1 to enable the full BAM count."
        fi
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

printf 'run\texit_status\treal_seconds\tdocker_client_user_seconds\tdocker_client_sys_seconds\tstdout_sha256\n' \
    > "$summary"

reference_stdout="$output_dir/reference.stdout"
cp "$reference_output" "$reference_stdout"
reference_checksum=$(checksum_file "$reference_stdout")
printf 'reference\t0\tNA\tNA\tNA\t%s\n' "$reference_checksum" >> "$summary"
echo "Using existing stock output: $reference_output"

stop_on_mismatch=${TELSEQ_STOP_ON_MISMATCH:-1}
had_failure=0

for thread_count in "${thread_grid[@]}"; do
    label="threads-${thread_count}"
    echo "Running $image with -t $thread_count..."
    if ! run_timed "$label" "${docker_run[@]}" \
        "${telseq_args[@]}" -t "$thread_count" "$container_bam"; then
        echo "FAIL: -t $thread_count exited unsuccessfully; see ${label}.stderr" >&2
        had_failure=1
        if [ "$stop_on_mismatch" != "0" ]; then
            break
        fi
        continue
    fi

    comparison_diff="$output_dir/${label}.diff"
    if comparison_mode=$("$output_comparator" \
        "$reference_stdout" "$output_dir/${label}.stdout" "$comparison_diff"); then
        if [ "$comparison_mode" = "exact" ]; then
            echo "PASS: -t $thread_count output is byte-identical to the reference"
        else
            echo "PASS: -t $thread_count inherited fields are byte-identical; candidate appends integer K"
        fi
        rm -f "$comparison_diff"
    else
        echo "FAIL: -t $thread_count output differs from the reference" >&2
        echo "Diff: $comparison_diff" >&2
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
echo "Only real_seconds measures the containerized TelSeq wall time."
echo "The user/sys columns measure the Docker client, not CPU time inside the container."
echo "Results: $output_dir"

exit "$had_failure"
