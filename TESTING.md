# Testing and benchmarking TelSeq Parallel

Correctness comes before speed. A parallel run should produce byte-identical
stdout to stock TelSeq when both programs receive the same BAM and the same
analysis parameters. Only the thread count should differ.

In particular, keep `-r`, `-k`, `-z`, `-e`, `-m`, `-u`, and `-w` identical.
The read-length option `-r` changes the number of `TEL` columns as well as the
calculation, so a mismatch there makes the outputs incomparable.

Published measurements are indexed in [benchmarks/README.md](benchmarks/README.md).

## Reproducible synthetic compatibility test

The repository includes a small synthetic BAM generator and an end-to-end
compatibility test. The fixture covers indexed-parallel edge cases:

- 64 reference sequences, including 8 empty references;
- 1,120 mapped alignments;
- one unmapped alignment that still has a reference coordinate;
- 8 unmapped, no-coordinate alignments;
- reverse-strand, duplicate, secondary, and failed-QC flags;
- a coordinate-sorted header and generated BAI index; and
- a telomeric final physical record, exercising stock TelSeq's legacy
  final-record contribution.

Build TelSeq, then run:

```bash
scripts/test_parallel_synthetic.sh \
    /path/to/stock/telseq \
    src/Telseq/telseq
```

On macOS, the bundled stock executable can normally be used:

```bash
scripts/test_parallel_synthetic.sh \
    bin/mac/telseq \
    src/Telseq/telseq
```

On Linux, use `bin/ubuntu/telseq` if it is compatible with the host, or supply
another unmodified TelSeq executable.

The test compares stock TelSeq with the new binary at `-t 1`, `-t 2`, `-t 22`,
and `-t 44`. Every stdout file must be byte-identical. It also checks that:

- `-t 22` starts 21 indexed workers plus one HTSlib compatibility scanner;
- `-t 44` starts 43 indexed workers plus one HTSlib compatibility scanner;
- the compatibility scanner fetches only the 8 indexed no-coordinate records,
  reports zero complete sequential scans, and retains the legacy EOF record;
  and
- the output contains `Total=1130`, `Mapped=1120`, and `Duplicates=3`.

The fixture contains 1,129 physical records. The expected `Total` of 1,130
preserves the original TelSeq final-record behavior.

The script prints its temporary artifact directory and retains the generated
BAM, index, stdout, stderr, timings, checksums, and any output differences.

The fixture is designed for correctness, not performance. It is too small to
measure useful parallel speedup.

## Synthetic thread-scaling regression

The Linux scaling regression creates a larger but highly compressible fixture,
warms both worker layouts, runs each layout three times, and compares median
wall time:

```bash
scripts/test_parallel_scaling.sh \
    src/Telseq/telseq \
    src/Test/generate_parallel_fixture
```

The test requires at least two logical CPUs. It verifies that `-t 2` and
`-t 4` produce byte-identical stdout and requires the median `-t 4` wall time
to be lower than the median `-t 2` wall time. The workload size can be changed
without editing the script:

```bash
TELSEQ_SCALING_READS_PER_REFERENCE=8000 \
TELSEQ_SCALING_READ_LENGTH=1000 \
    scripts/test_parallel_scaling.sh \
        src/Telseq/telseq \
        src/Test/generate_parallel_fixture
```

This is a regression test for working concurrency, not a substitute for a
representative WGS benchmark.

## Compare and benchmark on a real WGS BAM

The benchmark wrapper runs a grid of thread counts, captures timing and logs,
and compares every stdout file byte-for-byte with a stock reference.

### Run stock TelSeq as the reference

```bash
scripts/compare_and_benchmark.sh \
    /path/to/stock/telseq \
    /path/to/telseq-parallel \
    /path/to/sample.bam \
    1 2 4 8 16 22 44
```

When no thread grid is supplied, the default is `1 2 4 8 16 24`.

### Reuse existing stock output

If stock TelSeq has already been run on the exact same BAM with the exact same
analysis parameters, avoid another complete BAM scan:

```bash
scripts/compare_and_benchmark.sh \
    --reference-output /path/to/stock-result.tsv \
    /path/to/telseq-parallel \
    /path/to/sample.bam \
    22 44 80
```

The reference file must contain stock TelSeq stdout only. Do not use a file in
which stderr progress messages were mixed into the table.

### Forward TelSeq analysis parameters

Place analysis parameters after `--`. They are forwarded unchanged to every
parallel run:

```bash
scripts/compare_and_benchmark.sh \
    --reference-output /path/to/stock-result.tsv \
    /path/to/telseq-parallel \
    /path/to/sample.bam \
    22 44 80 \
    -- -k 7 -r 151 -u
```

The stock reference must have been produced using the same `-k 7 -r 151 -u`
arguments. Do not pass the following options after `--`:

- `-t` / `--threads`, because the wrapper controls the thread grid;
- `-f` / `--bamlist`, because the wrapper controls the input BAM; or
- `-o` / `--output-dir`, because the wrapper captures stdout.

### Benchmark configuration

The following environment variables alter wrapper behavior:

| Variable | Default | Effect |
| --- | --- | --- |
| `TELSEQ_BENCH_OUT` | timestamped directory | Select the artifact directory. |
| `TELSEQ_STOP_ON_MISMATCH` | `1` | Stop at the first failed run or output mismatch. Set to `0` to finish the grid. |
| `TELSEQ_RUN_BAM_COUNT` | `0` | Run the additional `samtools view -c` full BAM scan when set to `1`. |

For example:

```bash
TELSEQ_BENCH_OUT=benchmark-sample1 \
TELSEQ_STOP_ON_MISMATCH=0 \
    scripts/compare_and_benchmark.sh \
        --reference-output stock-sample1.tsv \
        src/Telseq/telseq \
        sample1.bam \
        1 8 16 22 44 \
        -- -r 151
```

The expensive `samtools view -c` count is deliberately disabled by default.
It adds a separate complete BAM pass and can take several minutes on WGS data.
Enable it only when the independent physical-record count is needed.

### Benchmark artifacts

Each run writes a timestamped directory containing:

| Artifact | Contents |
| --- | --- |
| `environment.txt` | Host, kernel, input paths, thread grid, TelSeq arguments, version information, and lightweight BAM checks. |
| `summary.tsv` | Exit status, real/user/system time, and stdout SHA-256 for every run. |
| `reference.stdout` or `stock.stdout` | Reference table used for comparisons. |
| `threads-N.stdout` | TelSeq result for thread count `N`. |
| `threads-N.stderr` | Progress log and `/usr/bin/time -p` output. |
| `threads-N.diff` | Unified diff created only when output differs. |

A successful comparison prints:

```text
PASS: -t 22 output is byte-identical to the reference
```

Matching `Total`, `Mapped`, and `Duplicates` columns alone is not sufficient.
Compare the complete stdout because `LENGTH_ESTIMATE`, all `TEL` bins, and all
`GC` bins participate in the result.

## How to benchmark fairly

Runtime depends heavily on storage and cache state. Record enough context to
make comparisons meaningful:

1. Use the same BAM, index, node, requested resources, and TelSeq parameters.
2. Confirm output identity before using a timing result.
3. Record whether the BAM was read from a warm operating-system cache or cold
   storage. Do not compare a cold first run with cached later runs as if thread
   count were the only difference.
4. Repeat useful thread counts when time permits and compare medians rather
   than a single run.
5. Watch CPU utilization, I/O throughput, and memory on the cluster. More
   threads can be slower once storage or decompression is saturated.
6. Include `-t 1` as the new binary's sequential baseline when possible.

The requested thread count is not a chromosome count. For `-t N`, one thread
is reserved for the HTSlib compatibility scanner and up to `N-1` workers
dynamically consume whole-reference tasks. More threads than useful reference
tasks do not create more parallel work.

The compatibility scanner uses the BAI to retrieve the no-coordinate tail
directly. Its progress log reports how many indexed records it fetched and
must report `full sequential scans: 0`. This retains no-coordinate alignments
and stock TelSeq's legacy final-record contribution without reading the whole
BAM a second time.

## Container validation

The release workflow builds the Linux AMD64 container from the current source.
As part of the Docker build it runs the full stock-compatibility fixture,
asserts indexed no-coordinate access, and runs the synthetic thread-scaling
regression. The image cannot be published unless all of these checks pass.

This container check complements rather than replaces comparison with an
unmodified stock TelSeq binary. The full synthetic script and at least one
representative real WGS comparison should be used before adopting new changes
in production.
