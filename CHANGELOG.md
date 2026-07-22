# Changelog

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
