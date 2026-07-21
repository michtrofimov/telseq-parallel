# Changelog

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
