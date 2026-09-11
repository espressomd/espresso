# P3M OpenMP parallel efficiency — analysis (2026-07-16)

Machine: Intel i7-12700K, **hybrid** 8 P-cores (@4.9 GHz) + 4 E-cores (@3.6 GHz),
20 hardware threads. Build: build-native (`-march=native`), current head.
10 000 particles, single MPI rank, OpenMP threads. Benchmark timed loop only.

## Scaling (1 → 4 OpenMP threads)

| config | t1 (ms) | t4 (ms) | speedup | efficiency | Amdahl serial fraction |
|---|---:|---:|---:|---:|---:|
| lj 10k              | 2.33 | 0.865 | 2.69× | 67% | ~16% |
| p3m 10k mesh=64     | 14.0 | 8.62  | 1.63× | 41% | ~48% |
| p3m 10k mesh=96     | 35.5 | 25.9  | 1.37× | 34% | ~64% |

Finer mesh scales *worse* — the added work (t1 +21.5 ms, t4 +17.3 ms for 64→96)
scales at only ~1.24×. So the mesh/FFT (reciprocal-space) path is what drags
P3M down; the particle-space path scales like LJ.

## Root cause: the reciprocal-space path does not use OpenMP threads

`perf stat`, p3m mesh=96:

| | t1 | t4 |
|---|---:|---:|
| CPUs utilized (task-clock / wall) | 1.00 | **1.80** |
| insn/cycle (P-core) | 2.07 | 1.73 |
| LLC-load-miss rate | 23% | 29% |

At 4 threads only **1.8 of 4 cores** are busy on average — threads sit idle
more than half the time. That is a serial-fraction signature, not a
bandwidth-only one (lj t4 reaches 2.8 CPUs on the same machine).

Flat `perf` profile of p3m mesh=96 at 4 threads is dominated by the FFT /
reshape / mesh-buffer work, not the particle kernels:

- `heffte::transpose_packer::unpack` + `fftw_cpy1d` + `fftw_cpy2d_pair`
  (the 3D-FFT stage transposes and copies)
- `__memmove` / `__memset` of the FFT/mesh buffers (8–13% each)
- `HostIterateTile<MDRangePolicy<Rank<3>>>` (the mesh-wise operations over the
  96³ grid — influence function, etc.)
- particle kernels (`SpecializedForcesKernel`, `link_cell`, `init_forces`)
  together only ~10%

ESPResSo parallelizes P3M's mesh/FFT over **MPI ranks** (the mesh is
domain-decomposed and heFFTe does distributed transposes). Within a single
rank the FFT (serial FFTW plans) and the heFFTe transposes get essentially no
thread parallelism, and the mesh element-wise ops that *do* thread
(`MDRangePolicy`) are memory-bandwidth-bound. So OpenMP threads speed up only
the particle work; the reciprocal-space half (which grows with mesh size)
stays serial. Hence 64% serial at mesh=96.

## Amplifying factors

1. **Hybrid cores + no pinning — MINOR.** Pinning the 4 threads to 4 distinct
   P-cores (`taskset -c 0,2,4,6`) gave p3m mesh=96 t4 = 25.4 ms vs 25.9 ms
   unpinned (lj: 0.888 vs 0.865) — within noise. So the E-core/HT imbalance is
   not the cause; the serial fraction dominates wherever the threads run.
   Decisive check: **8 threads on 8 P-cores = 26.2 ms, no better than 4** — more
   threads do nothing because the reciprocal-space path is serial (Amdahl
   ceiling at s≈0.64 is ~1.5× for any thread count).

2. **Memory bandwidth.** The mesh ops and FFT transposes are bandwidth-bound
   (LLC-miss 23%→29%, IPC 2.07→1.73 from t1→t4), so even the threaded mesh work
   scales sublinearly.

3. **Fork/join over many small regions.** Even lj is only 67% efficient (~16%
   serial): each step runs many small Kokkos `parallel_for` regions (force
   init, verlet build, pair loop, ScatterView reduce, two integrator steps),
   and at 10k particles/thread the per-region fork/join + ScatterView duplicate
   reduction is a fixed overhead that does not shrink with threads.

## Threaded FFTW — tested, does NOT help

heFFTe 2.4.1 has no OpenMP: no `#pragma omp` in its reshape/pack code, and its
FFTW executor creates plans with `fftw_plan_many_dft` and never calls
`fftw_plan_with_nthreads` (that call exists only in heFFTe's own `speed3d`
benchmark, not the library). The build DOES link `libfftw3_omp` (HeffteConfig),
so the threaded planner is available but never switched on.

Experiment: an `LD_PRELOAD` shim called `fftw_init_threads()` +
`fftw_plan_with_nthreads(4)` at startup (verified engaged:
`fftw_planner_nthreads() == 4`), so heFFTe's plans were built with the 4-thread
planner. Result — **no change**:

| config | FFTW 1 thread | FFTW 4 threads |
|---|---:|---:|
| p3m mesh=96, 4 Kokkos thr | 25.6 ms | 25.8 ms |
| p3m mesh=64, 4 Kokkos thr | 8.38 ms | 8.62 ms |

Why: the single-rank serial cost is dominated by heFFTe's **reshape/transpose
data movement** (`transpose_packer::unpack`, `__memmove`, `__memset` of the mesh
buffers), which is heFFTe serial code with no OpenMP — not the FFTW transforms.
Threaded FFTW only touches the FFT-execute (a minority), and the batched
size-96/128 1D FFTs are too small for FFTW to thread profitably anyway.

## Prototype: OpenMP in heFFTe's reshape kernels — WORKS

Since threaded FFTW addressed the wrong part, I prototyped OpenMP on heFFTe's
own reshape/pack kernels (`heffte_pack3d.h`): `#pragma omp parallel for` on the
outer loop of `direct_packer::pack`/`unpack`, `collapse(2)`/`collapse(3)` on the
six `transpose_packer::unpack` branches, and the post-FFT `data_scaling::apply`.
The `pack` inner iterator was rewritten to an explicit destination offset so the
outer loop parallelizes. Each output element is written exactly once from a
fixed source, so the transform result is a pure permutation — order-independent,
hence bitwise-identical. Patch: `.superpowers/sdd/heffte-reshape-openmp.patch`
(applied to build-native's FetchContent heFFTe copy for the prototype).

Result (build-native, 10k particles):

| config | t1 | t4 before | t4 after | speedup | efficiency |
|---|---:|---:|---:|---:|---:|
| p3m mesh=96 | 36.1 ms | 25.9 ms | **22.5 ms** | 1.37×→**1.61×** | 34%→**40%** |
| p3m mesh=64 | 14.0 ms | 8.62 ms | **7.14 ms** | 1.62×→**1.96×** | 41%→**49%** |

4-thread P3M is 13–17% faster; 1-thread unchanged. Native identity bitwise
preserved (lj `1671c333`, p3m `d2dfa428`) — confirmed the permutation argument.
The remaining serial fraction is the single-threaded FFTW execute plus the
reshape3d buffer management (memmove/memset) that this patch does not touch, so
efficiency is up but not near 100%.

Integration path (not done — heFFTe is a FetchContent dependency): either apply
the patch via a `PATCH_COMMAND` in ESPResSo's heFFTe FetchContent block, or
upstream the pragmas to heFFTe (the cleaner option — heFFTe would want them
guarded/uniform across backends).

## What would actually help (directions)

- **Use MPI ranks, not OpenMP threads, for P3M** — heFFTe's mesh/FFT is built to
  scale over ranks (distributed transposes). The single-rank multi-thread mode
  is the wrong axis for the reciprocal space; threads are for the particle work.
- **Threading heFFTe's reshapes** would need patching/upstreaming heFFTe itself
  (its transpose/pack kernels are serial); threaded FFTW alone is not enough
  (shown above).
- Thread pinning: marginal here (shown above).
- Particle side: fewer/larger fused Kokkos regions + a lighter force reduction
  would lift the ~67% lj ceiling, but that is a smaller prize than P3M's
  reciprocal-space serial fraction.

## Addressing the remaining gap (2 more commits)

Profiling the packer-threaded build (p3m mesh=96, 4 thr) showed the next serial
slices were `__memmove`/`__memset` of mesh-sized buffers. A 1/2/4-thread curve
(mesh96: 1->2 1.32x, 2->4 1.27x) showed residual threadable headroom, not a hard
bandwidth wall.

- `a77bee9c6d` heFFTe reshape OpenMP (via cmake/heffte.patch PATCH_COMMAND):
  mesh96 t4 25.9->21.1 ms.
- `dfa71afb60` p3m mesh crop/pad copies (extract_block, pad_with_zeros_discard_imag,
  in-tree): explicit-offset outer loop under `omp parallel for`; mesh96 t4
  21.1->19.9 ms, mesh64 t4 7.3->6.7 ms. Bitwise-identical (permutation).

Cumulative single-rank P3M @ 10k/4 thr: mesh=96 25.9->19.9 ms (efficiency
34%->45%), mesh=64 8.62->6.7 ms (41%->51%); 1-thread unchanged; identity
bitwise-preserved throughout.

Still remaining (diminishing returns, not pursued): the std::vector zero-init
memset in crop/pad (extract_block's is fully overwritten -> could allocate
uninitialized), the single-threaded FFTW execute (threading it gave nothing),
and the bandwidth-bound mesh MDRange ops. The real lever beyond this is MPI
ranks (heFFTe's design distributes the mesh + memory bandwidth across ranks).

## MPI-rank strong scaling (2026-07-17) — the intended axis

Confirmation that MPI ranks, not OpenMP threads, are the lever for P3M.
Fixed total N = 16 000, ranks bound one-per-P-core (`mpirun --bind-to core
--map-by core`), `OMP_NUM_THREADS=1`, mesh fixed via `--mesh`. Machine quiet
(1-3% foreign) for all 12 runs. build-heffte-omp at head.

| config | t1 | t2 (eff) | t4 (eff) | t8 (eff) |
|---|---:|---:|---:|---:|
| LJ 16k          | 3.90 ms | 2.31 (84%) | 1.45 (67%) | 0.94 (52%) |
| P3M 16k mesh=64 | 20.4 ms | 11.95 (86%)| 6.96 (74%) | 5.01 (51%) |
| P3M 16k mesh=96 | 39.9 ms | 26.7 (75%) | 14.7 (68%) | 12.5 (40%) |

**Axis comparison, P3M mesh=96, 4-way parallel:** OpenMP 4 threads = 45.8%
efficiency (20.0 ms/step); MPI 4 ranks = **67.9%** (14.7 ms/step). MPI ranks
deliver the parallelism OpenMP could not, because heFFTe domain-decomposes the
mesh + distributed transposes across ranks — exactly the reciprocal-space work
that stayed serial under threads. 8 ranks reach 12.5 ms/step (3.19x over 1
rank), well past the best single-rank threaded result.

**Where rank scaling tapers:** mesh=96 efficiency falls to 40% at R=8 — with
only 2000 particles/rank and a 96^3 mesh split 2x2x2, heFFTe's all-to-all
transpose becomes communication-bound. mesh=64 holds ~74% through R=4. So on
this single node the sweet spot is ~4 ranks; more ranks help absolute latency
(12.5 vs 14.7 ms) but at falling efficiency. Real clusters push the crossover
out with faster interconnect + more mesh per rank.

Conclusion: the earlier single-rank OpenMP work (heFFTe packer + crop/pad
threading) was the right micro-optimization for the 1-rank case, but the
order-of-magnitude lever is MPI ranks, which is heFFTe's design point.
