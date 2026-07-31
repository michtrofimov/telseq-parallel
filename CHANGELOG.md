# Changelog

## Unreleased

## 0.5.4 - 2026-07-31

### Changed

- the appended `-u` row now reports position-matched pipe-delimited read
  groups, libraries, and samples, such as `RG1|RG2`, `LIB1|LIB2`, and
  `SAMPLE1|SAMPLE2`;
- reads with missing `RG` tags or tags not declared in the BAM header are
  skipped and do not contribute to either regular or appended `-u` results.

## 0.5.3 - 2026-07-31

### Changed

- `-u` now retains every declared read-group row and appends one exact
  whole-BAM aggregate row, matching the additive output style of `-m`;
- the appended `-u` row explicitly reports `UNKNOWN` for read group, library,
  and sample, while the regular rows preserve their BAM-header metadata;
- reads with missing or undeclared `RG` tags still contribute to the `-u`
  whole-BAM aggregate.

## 0.5.2 - 2026-07-29

### Changed

- restored the `v0.5.0` read-group behavior: default output contains regular
  rows only, while `-m` retains those rows and appends one merged
  weighted-average row.

## 0.5.1 - 2026-07-28

### Changed

- BAMs with more than one read group append the merged weighted-average row by
  default while retaining every regular read-group row;
- `-m` remains accepted for command-line compatibility and produces the same
  output as the new default behavior.

## 0.5.0 - 2026-07-28

### Changed

- result tables, including `-m` merged read-group rows, are written to stdout
  as well as to the selected `--output-tsv` file;
- `-m` retains the regular per-read-group rows and appends the merged
  weighted-average row instead of replacing those regular rows;
- `--output-tsv /dev/stdout` still emits exactly one clean TSV copy and keeps
  progress logging off that stream.

## 0.4.1 - 2026-07-28

### Changed

- file-output runs stream the complete progress and reference-profile log to
  standard output in real time while also saving it to `--output-log`;
- direct result streaming with `--output-tsv /dev/stdout` keeps stdout
  byte-clean and continues routing diagnostics only to the selected log.

## 0.4.0 - 2026-07-28

### Changed

- file-output runs report concise start and completion messages, including the
  result TSV and profile-log paths, to standard output; direct TSV streaming to
  `/dev/stdout` remains byte-clean.

## 0.3.3 - 2026-07-28

### Added

- `--output-tsv` as a descriptive alias for the inherited `-o` and
  `--output-dir` result-file option;
- `--output-log` for writing progress, diagnostics, timing, and mapped-reference
  profile records to a dedicated log file;
- single-BAM runs derive `<BAM basename>.telseq.tsv` and
  `<BAM basename>.telseq.log` in the current working directory when paths are
  not supplied.

### Changed

- mapped-reference profiling is enabled by default for every parallel run;
- multi-BAM, bamlist, and piped-input runs require explicit TSV and log paths
  because they do not have one unambiguous BAM basename;
- benchmark and synthetic-test helpers explicitly route the new artifact
  destinations through standard output and error when capture is required.

### Compatibility

- `-o`, `--output-dir`, and the explicit `--profile-references` flag remain
  accepted;
- one-BAM runs now create result and log artifacts instead of using standard
  output and error as their primary destinations.

## 0.3.2 - 2026-07-23

### Added

- append an integer `K` column to every result row so the effective automatic
  or explicit telomeric-repeat threshold is recorded with the estimate.
- stock-output comparison helpers validate and remove only the additional
  final `K` column before requiring every inherited output byte to match.

### Fixed

- require explicit `-k` values to be complete integers instead of accepting
  an integer prefix from values such as `10.5` or `10x`.

### Compatibility

- result calculations and inherited columns remain unchanged, but raw stdout
  is no longer byte-identical to unmodified stock TelSeq because of the
  additional final `K` column.

## 0.3.1 - 2026-07-23

### Changed

- when `-k` is omitted, select the smallest motif-repeat count covering at
  least 40% of the configured read length; explicit `-k` values are unchanged.

### Compatibility

- the automatic threshold remains `k=7` for the default 100-base reads, but
  differs from stock TelSeq's fixed implicit `k=7` at longer read lengths;
  pass an explicit `-k 7` when reproducing that legacy behavior.

### Documentation

- publish the real-WGS version 0.3 reference-window scaling benchmark, where
  wall time reached 88.28 seconds at 80 requested threads and all measured
  thread counts produced the same result hash.

## 0.3.0 - 2026-07-22

### Added

- `compare_and_benchmark_docker.sh` for timing released container images and
  comparing every thread count byte-for-byte against saved stock output.
- Docker benchmarks preserve the host BAM path inside the container so
  inherited path text in stdout does not cause a false comparison failure.
- `--profile-references` emits per-reference scheduler assignment, read counts,
  and timing offsets to stderr for diagnosing parallel scaling limits without
  changing result stdout.
- successful `master` builds publish a moving GHCR development image for
  cluster profiling without creating a numbered release.
- long references are divided into dynamically scheduled indexed windows,
  with alignment-start ownership preventing duplicate counts for reads that
  overlap adjacent windows;
- `--reference-window-size` tunes window granularity or disables splitting for
  direct comparison with the previous whole-reference scheduler;
- synthetic boundary coverage verifies byte-identical output for reads that
  end at, span, start at, or start immediately after a window boundary.
- mapped-reference tasks are prioritized by BAI record-count estimates so
  dense short references start early instead of becoming late stragglers;
- reference profiles include the record-count estimate used by the scheduler.
- `--primary-chromosomes-only` restricts analysis to exact human autosomes
  1-22 and sex chromosomes X/Y (with optional `chr` prefix), excluding all
  contigs, mitochondrial references, and no-coordinate reads;
- strict primary-chromosome mode disables the compatibility scanner and makes
  every requested thread available to indexed primary-reference windows;
- synthetic tests cover unprefixed and `chr`-prefixed primary names, reject
  alt/mitochondrial/contig distractors, and compare serial and parallel output.

## 0.2.0 - 2026-07-22

### Added

- HTSlib-based direct indexed retrieval of no-coordinate BAM records;
- scanner instrumentation reporting indexed records fetched and confirming
  that no complete sequential scan was performed;
- a release-gating test that proves the compatibility scanner reads only the
  synthetic fixture's no-coordinate tail;
- a synthetic scaling regression comparing median wall time at `-t 2` and
  `-t 4` while requiring byte-identical output.

### Changed

- the compatibility worker no longer reads the entire BAM alongside the
  indexed reference workers;
- when a BAM has no no-coordinate tail, only its highest populated reference
  is inspected to recover stock TelSeq's legacy final-record contribution;
- parallel mode now requires a standard BAI readable by both BamTools and
  HTSlib, including no-coordinate-tail metadata normally written by
  `samtools index`;
- the Linux AMD64 container includes the HTSlib runtime.

### Compatibility

- the counting and tabular output calculation remain unchanged;
- stock TelSeq's final-record contribution is intentionally retained;
- `-t 1` continues to use the inherited sequential implementation.

## 0.1.0 - 2026-07-21

First parallel TelSeq release.

### Added

- `-t` / `--threads` for indexed parallel scanning of one coordinate-sorted
  BAM;
- dynamic whole-reference scheduling across independent BamTools readers;
- a compatibility scanner for no-coordinate alignments and stock TelSeq's
  legacy final-record contribution;
- byte-for-byte comparison and timing against existing stock output;
- forwarding of TelSeq analysis arguments through the benchmark wrapper;
- a reproducible 64-reference synthetic BAM compatibility test;
- a Linux AMD64 container image published to GHCR.

### Compatibility

- `-t 1` retains the original sequential scan path;
- parallel output is designed to preserve stock TelSeq counting behavior;
- `-t > 1` requires a coordinate-sorted BAM and a readable BAI/BTI index.

### Known limitation

The compatibility scanner performs one complete sequential BAM pass. Indexed
workers collectively perform another pass over reference-assigned records, so
storage bandwidth and decompression can limit scaling at high thread counts.
