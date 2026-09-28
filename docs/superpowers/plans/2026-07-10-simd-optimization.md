# SIMD Optimization Campaign (march=native)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Three independent SIMD workstreams on the component-major ParticleStore, evaluated under `-march=native` (AVX2+FMA on this host), each kept ONLY if it meets the acceptance bar.

**Acceptance bar (user-set, per workstream):** ≥3% improvement on ≥2 of the six benchmark configurations, OR ≥5% on one configuration with no slowdown (>1% regression = slowdown) on any other. Measured by order-balanced same-flags A/B (≥6 reps, alternating first runner) against the pre-workstream HEAD built with identical `-march=native` flags. perf (hardware sampling) guides the work; the A/B decides keep/discard.

**Correctness per workstream:** unit tests ALL; the physics pattern battery; same-flags bitwise identity where the rewrite preserves per-element arithmetic order — where batching legitimately reorders arithmetic (reductions), document it and gate via the physics/statistical tests instead. Standing compile gates (ci-mirror, gcc14, maxset, tidy, dox) stay on the CANONICAL flags (no march=native) — the campaign must not break them.

**Build setup:** `-march=native` is evaluated in dedicated build dirs (`/ssd/weeber/eliminate_Particle_builds/build-native` for the working tree; per-round reference builds from pinned worktrees). The canonical `build/` config (and the recorded baselines/identity reference) stay unchanged; whether `-march=native` becomes a recommended configuration is decided at the end from the data.

## Workstreams

### WS1: Batched minimum-image (non-Lees-Edwards fast path)
- `get_mi_vector_batch`: given one reference position and N neighbor positions (SoA spans), produce N minimum-image vectors — branchless per-axis periodic fold, vectorizable (no Lees-Edwards: guard with the existing LE-inactive predicate and fall back to the scalar path when LE is active).
- Use in the Cabana pair kernel (ForcesKernel): batch the per-neighbor MI vectors for each particle's Verlet neighbors.
- `get_mi_dist2_batch`: N squared minimum-image distances; used in the linked-cell Verlet BUILD as a cheap vectorized gate over each cell-pair's candidates BEFORE the scalar verlet criterion runs on survivors.
- Arithmetic-order note: per-pair values identical (no reduction) → same-flags bitwise identity expected for the build gate; the pair-kernel force ACCUMULATION order must be preserved (batch the MI computation, keep the existing accumulation loop order).

### WS2: SIMD-friendly P3M charge/force assignment
- Rewrite charge assignment (particle→mesh) and force interpolation (mesh→particle) with SIMD-friendly inner loops (contiguous mesh-line accesses, hoisted weights); consider `Kokkos::ScatterView` for the charge-assignment accumulation instead of the custom reduce.
- Reduction-order WILL change if ScatterView is adopted → p3m identity legitimately changes; gate via coulomb/p3m test battery (energies/forces vs tolerances) + document.

### WS3: SIMD-friendly force init + ghost-force init incl. Langevin
- Restructure `init_forces_and_thermostat` (and the ghost force reset) as column kernels with contiguous component loops (component-major columns make x/y/z runs contiguous); Langevin friction+noise vectorized where the Philox draws allow (per-particle id-keyed draws stay scalar; the deterministic friction arithmetic vectorizes).
- Per-particle values unchanged, no reductions → same-flags bitwise identity expected.

## Protocol per workstream
1. perf-profile the target slot under march=native (steady-state, -D skip).
2. Implement behind the smallest reasonable seam (no config toggle unless the LE guard requires one).
3. Correctness gates; then order-balanced A/B all six configs vs the pre-workstream reference binary (same flags).
4. KEEP if the bar is met; else REVERT (keep findings in the report). Ledger each verdict.
5. Canonical-flag compile gates + identity (canonical build unaffected) before any push. Push only to `me`.
