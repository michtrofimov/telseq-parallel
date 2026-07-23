# Changelog

## Unreleased

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
