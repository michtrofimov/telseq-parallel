# Real WGS v0.3 reference-window scaling benchmark — 2026-07-22

This benchmark measured the version 0.3.0 Docker image on one real WGS BAM at
4, 8, 12, 23, 46, and 80 requested threads. It used default compatibility mode,
not `--primary-chromosomes-only`, and therefore measures the reference-window
scheduler while retaining stock-compatible whole-BAM counting. The same BAM
and infrastructure were used for the version 0.2 benchmark. The sample
identifier and infrastructure details are omitted.

![TelSeq Parallel v0.3 wall time versus requested threads](wall-time-vs-threads.svg)

## Results

| Threads | Mapped workers¹ | Wall time (s) | Wall time (min) | Speedup vs. `-t 4` | Reduction |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 3 | 1296.86 | 21.61 | 1.00× | 0.00% |
| 8 | 7 | 576.65 | 9.61 | 2.25× | 55.53% |
| 12 | 11 | 367.72 | 6.13 | 3.53× | 71.65% |
| 23 | 22 | 190.84 | 3.18 | 6.80× | 85.28% |
| 46 | 45 | 106.80 | 1.78 | 12.14× | 91.76% |
| 80 | 79 | **88.28** | **1.47** | **14.69×** | **93.19%** |

¹ In default mode, one requested thread is the HTSlib compatibility scanner
and up to `N-1` threads process dynamically scheduled mapped-reference windows.

The raw measurements are in [summary.tsv](summary.tsv).

## Interpretation

Wall time fell from 21.61 minutes at `-t 4` to 1.47 minutes at `-t 80`, a
14.69× speedup and a 93.19% reduction within this measured series. Scaling
continues beyond 23 threads: increasing the request from 23 to 80 threads
saved 102.56 seconds, or 53.74%.

This is the expected effect of the version 0.3 scheduler. Long references are
split into indexed windows and placed in the same dynamic work queue as short
references. Workers no longer wait for one thread to finish an entire long
chromosome, while alignment-start ownership prevents boundary reads from
being counted twice. BAI record-count estimates are used to start expensive
windows earlier and reduce late stragglers.

### Comparison with version 0.2

| Threads | v0.2 whole-reference (s) | v0.3 windows (s) | v0.3 reduction |
| ---: | ---: | ---: | ---: |
| 4 | 1309.12 | 1296.86 | 0.94% |
| 8 | 596.09 | 576.65 | 3.26% |
| 12 | 392.53 | 367.72 | 6.32% |
| 23 | 342.62 | 190.84 | 44.30% |
| 46 | 337.16 | 106.80 | 68.32% |
| 80 | 335.88 | 88.28 | 73.72% |

The small difference at low thread counts and large difference at high thread
counts are consistent with removing the whole-chromosome scheduling ceiling.
The earlier version 0.2 curve plateaued near 23 requested threads; the v0.3
curve continues improving through the largest tested setting.

## Output comparison

Every container run exited successfully and produced the same SHA-256:

```text
23a1123550574145e5d3e60fa706c350c576894d13b4e432f374148f1a78c969
```

This is also the result hash produced by every parallel run in the version 0.2
benchmark. The saved reference file has a different whole-file hash because
its header differs; the result rows are unaffected.

## Limitations

- There is one timing observation per thread count.
- No fresh `-t 1` measurement was included, so speedups use `-t 4` as the
  baseline.
- Cache state, node model, filesystem throughput, container CPU utilization,
  and storage load were not recorded with the submitted summary.
- This series tests default whole-BAM compatibility mode. It does not measure
  the optional `--primary-chromosomes-only` filter.
- Results should not be generalized to a different BAM or storage system
  without benchmarking it directly.
