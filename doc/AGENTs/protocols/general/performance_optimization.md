---
description: Optimize GPU kernels, memory layouts, and numerical precision while preserving correctness via regression-guarded benchmarking.
---

# Performance Optimization Protocol

## Purpose

Maximize throughput (steps/second, evaluations/second) of a compute kernel or solver subject to a fixed correctness constraint. The loop tries parameter variants, benchmarks each, and keeps the fastest that does not break validation.

This protocol applies to:
- GPU kernel tuning (workgroup size, memory layout, vector width)
- Algorithmic variants (atomics vs. reduction, direct vs. grid-based)
- Precision tradeoffs (double → single → half)

---

## Agentic Loop Integration

```
qualitative validation passes
    |
    v
THIS PROTOCOL (Level 6)
    |
    v
human review (if speedup targets not met)
```

**Never optimize broken code.** The full correctness suite (compile → topology → parity → qualitative) must pass before any timing measurement is trusted.

---

## Correctness Regression Guard

For every candidate variant, re-run:

1. **Level 0–1**: Compile, no crash, finite outputs.
2. **Level 2**: Topology unchanged (same neighbor lists, same padding).
3. **Level 3–4**: Parity against the *previous best* variant (not the original reference). Tolerance: `method_class` strategy.

Only if all pass, record the timing. If any fail, **discard the variant immediately** and log the failure mode.

---

## Benchmarking Procedure

### 1. Warmup & Statistical Timing

```
for each variant:
    # Warmup (exclude from timing)
    for i in range(n_warmup):
        run(variant)

    # Timed runs
    times = []
    for i in range(n_runs):
        t0 = timer()
        run(variant)
        queue.finish()   # GPU: ensure kernel completed
        times.append(timer() - t0)

    report median, min, max, stddev
```

Use **median**, not mean, to suppress outlier jitter. `n_warmup = 10`, `n_runs >= 30`.

---

### 2. Metrics to Report

| Metric | How to Measure | Target |
|--------|---------------|--------|
| Kernel time | OpenCL event profiling / CUDA events | Minimize |
| End-to-end throughput | `N_evaluations / total_time` | Maximize |
| Memory bandwidth | Bytes read+written / kernel_time | Compare to theoretical peak |
| Occupancy | Active warps / max warps | > 50% |
| Host-device transfer | `queue.finish()` around upload/compute/download | Minimize |
| Energy efficiency | throughput / watts (if available) | Maximize |

**Roofline model**: Plot achieved bandwidth vs. arithmetic intensity to determine if the kernel is memory-bound or compute-bound. This guides which optimizations to try next.

---

### 3. Parameter Sweep Space

| Parameter | Values to Try | Impact |
|-----------|--------------|--------|
| Workgroup size | 32, 64, 128, 256, 512 | Occupancy, register pressure, barrier cost |
| Memory layout | AoS → SoA | Coalescing, cache efficiency |
| Vector width | scalar, float2, float4, float8 | Instruction throughput |
| Precision | double, single, mixed | Speed vs. accuracy |
| Reduction strategy | atomic_add, local reduction, recoil assembly | Contention vs. extra passes |
| Loop unrolling | 1×, 2×, 4×, 8× | Instruction cache, register pressure |
| Local memory caching | none, neighbor positions, parameter tables | Reuse vs. capacity |

**Sweep strategy**: Grid search for coarse bounds, then local search around best. Do not try all combinations — use physical reasoning to prune (e.g., if memory-bound, try SoA first; if compute-bound, try unrolling first).

---

### 4. Decision Criteria

A variant is accepted as the new best if:

```
median_time < 0.95 * best_time          # At least 5% faster
AND
parity_pass(variant, best_variant)      # Correctness preserved
```

If margin is < 5%, keep the simpler variant (fewer lines, fewer branches, more readable).

---

## Performance Artifact

```
debug/<date>_<task>/
  ├── benchmark_summary.txt
  │       Variant        Median(ms)  Min    Max    Std    Occupancy  Parity
  │       baseline       12.34       11.9   13.2   0.4    0.72       —
  │       wg_128_soa     8.91        8.7    9.4    0.2    0.85       PASS
  │       wg_256_unroll  9.12        8.9    9.6    0.3    0.68       PASS
  │       wg_64_atomics  15.23       14.8   16.1   0.5    0.91       PASS
  │       → Best: wg_128_soa (1.38x speedup)
  ├── roofline.png
  └── best_variant.diff   # Code diff vs. baseline
```

---

## Rules

1. **Speed is worthless if correctness is lost.** The regression guard is absolute.
2. **Use event profiling, not wall-clock, for kernel time.** Wall-clock includes host overhead.
3. **Report median + min/max, not just average.** GPU timing has long-tail jitter.
4. **If a variant fails correctness, log why before discarding.** The failure mode guides the next attempt (e.g., "atomic collision on recoil → try local reduction").
5. **Never optimize for a single system size.** Test scaling: small (10 atoms), medium (100), large (1000+) to catch cache/occupancy cliffs.
