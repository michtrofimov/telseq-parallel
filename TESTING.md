# Parallel TelSeq testing

## Reproducible synthetic test

The repository includes a synthetic BAM generator and an end-to-end
compatibility test. The generated fixture is deliberately small, but covers
the indexed-parallel edge cases:

- 64 reference sequences, including 8 empty references;
- 1,120 mapped alignments;
- one unmapped alignment that still has a reference coordinate;
- 8 unmapped, no-coordinate alignments;
- reverse-strand, duplicate, secondary, and failed-QC flags;
- a coordinate-sorted header and a generated BAI index;
- a telomeric final physical record, exercising stock TelSeq's legacy
  final-record contribution.

After building TelSeq, run:

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
the path to another unmodified TelSeq executable.

The test generates the BAM and index in a temporary directory, then compares
stock TelSeq with `-t 1`, `-t 2`, `-t 22`, and `-t 44`. Every stdout file must
be byte-identical. It also verifies that:

- `-t 22` starts 21 indexed workers plus one compatibility scanner;
- `-t 44` starts 43 indexed workers plus one compatibility scanner;
- the output contains `Total=1130`, `Mapped=1120`, and `Duplicates=3`.

There are 1,129 physical records in the fixture. The expected total of 1,130
preserves stock TelSeq's legacy contribution of the final alignment.

The script prints the temporary artifact directory and keeps the generated
BAM, index, stdout, stderr, timings, and checksums for inspection.

## Real WGS compatibility and timing

If stock TelSeq output already exists for the exact same BAM path and options,
avoid another stock scan:

```bash
scripts/compare_and_benchmark.sh \
    --reference-output /path/to/stock-result.tsv \
    src/Telseq/telseq \
    /path/to/sample.bam \
    22 44 80 \
    -- -k 7 -r 100 -u
```

Arguments following `--` are forwarded unchanged to every new TelSeq run. The
reference output must have been produced using those same arguments. The
wrapper reserves `-t`, `-f`, and `-o` because it controls threads, BAM input,
and stdout capture itself.

The benchmark skips `samtools view -c` by default. It should be enabled with
`TELSEQ_RUN_BAM_COUNT=1` only when the additional full BAM pass is useful.

The synthetic fixture is intended to test correctness and concurrent task
layout. It is too small to measure useful speedup; performance must be measured
on a representative BAM and storage system.
