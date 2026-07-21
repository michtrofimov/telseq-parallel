# Real WGS thread-scaling benchmark — 2026-07-21

This benchmark measured one real WGS BAM at 1, 4, 8, 12, 22, 44, and 80
requested threads. The sample identifier and infrastructure details are omitted.
The 1–12 and 22–80 points came from two invocations of the benchmark wrapper.

![TelSeq Parallel wall time versus requested threads](wall-time-vs-threads.svg)

## Results

| Threads | Wall time (s) | Wall time (min) | Speedup vs. `-t 1` | Approx. CPU equivalents¹ |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3277.24 | 54.62 | 1.00× | 1.00 |
| 4 | **2288.93** | **38.15** | **1.43×** | 2.69 |
| 8 | 2307.01 | 38.45 | 1.42× | 2.71 |
| 12 | 2311.14 | 38.52 | 1.42× | 2.76 |
| 22 | 2313.10 | 38.55 | 1.42× | 2.77 |
| 44 | 2312.74 | 38.55 | 1.42× | 2.75 |
| 80 | 2306.38 | 38.44 | 1.42× | 2.75 |

¹ `(user_seconds + sys_seconds) / real_seconds`; this is an average over the
whole run, not a direct CPU-utilization measurement.

The raw measurements are in [summary.tsv](summary.tsv).

## Interpretation

`-t 4` was the fastest observed setting. It saved 988.31 seconds, or 16.47
minutes, relative to `-t 1`: a 30.16% wall-time reduction and 1.43× speedup.

All measurements from 4 through 80 threads fall within 24.17 seconds of one
another, a spread of only 1.06% relative to the fastest run. The requested
thread count increases twenty-fold across that range, while average consumed
CPU remains close to 2.7 CPU equivalents. On this system, additional threads
therefore increased resource demand without producing additional throughput.

This plateau is consistent with the full compatibility scan becoming the
critical path once three indexed workers are available, together with limits
from BAM decompression or storage throughput. The timing table alone cannot
distinguish those causes; I/O and CPU profiling would be needed to prove the
bottleneck.

For this BAM and infrastructure, `-t 4` is the best practical choice among the
tested settings because it has the lowest observed wall time and requests far
fewer CPUs. Repeated measurements of `-t 4` and a nearby setting should be used
to confirm the choice if small timing differences matter.

## Output comparison

Every TelSeq Parallel run exited successfully and produced the same SHA-256:

```text
cc4ea6e2046822960009b0cc3d73cb4310274afe35c0dde6fbffc4d2a8af9a12
```

The stored stock reference has a different whole-file SHA-256:

```text
17bdc6c5d872c397d35d31c8be42b7315f7ba399ebfb25b3157b367c51777626
```

Manual inspection reported identical data rows with a header-only difference.
Because the original files are not stored here, that normalized comparison has
not been independently reproduced from this repository. Future benchmark runs
should retain the diff or a header-normalized checksum alongside the summary.

## Limitations

- There is one timing observation per thread count.
- Cache state, BAM size, node model, filesystem, and storage load were not
  recorded with the submitted summary.
- The measurements came from two wrapper invocations.
- These results should not be generalized to a different BAM or storage system
  without benchmarking it directly.
