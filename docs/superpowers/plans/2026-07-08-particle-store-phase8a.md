# ParticleStore Migration — Phase 8a Implementation Plan (de-proxy CPU integrator/thermostat hot loops)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite the per-step integrator and thermostat hot loops as direct column kernels — column pointers hoisted once per loop, no per-row Particle-view rebinding, vectorizer-transparent bodies — closing the lj-1rank proxy-tax residual and putting the ORIGINAL 5% budget back in reach for the phase-8 checkpoint.

**Architecture:** First of four phase-8 sub-projects (8a → 8b LB coupling → 8c GPU single-rank device path → 8d cleanup; exploration `.superpowers/sdd/phase8-exploration.md` is the basis — read it before each task). Ground truth: the hot loops already run `Kokkos::parallel_for` over `[0,n_local)` via `for_each_local_particle` (`for_each_particle.hpp:37-62`), but each iteration rebinds a `Particle` view and every accessor pays `view_host()` + address recompute + stride, opaque to the vectorizer. De-proxying = same loop, same iteration order, same arithmetic per element, with hoisted raw column pointers — BITWISE identity must hold (no reductions in these loops; RNG is per-particle keyed — verify the key derivation is untouched). Scope: velocity-Verlet translation step1/step2, VV rotation, Langevin (translation + rotation), Brownian (same family, same treatment) — NpT/Stokesian/RATTLE stay on the view path (identity-safe, not hotspots, per exploration §6).

**Spec:** section 3 phase 8. Perf target: lj1 ≤1.05 vs phase-0 at the 8a gate (the amendment retires at the phase-8 checkpoint if all configs ≤5%); acceptance protocol = order-balanced min-of-means A/B vs the phase-7 tip binary (build once at `/ssd/weeber/eliminate_Particle_builds/ab-phase7`) + the sequential gate.

## Global Constraints

- All standing constraints (worktree; `make -j8`; controller benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`).
- MANDATORY before every commit: `maintainer/format/clang-format.sh` on changed files.
- Stale-ABI rule: header change → clean-rebuild `espresso_core` + `unit_tests_executables` everywhere tested; NFS → `cmake --fresh` absolute paths.
- CI gates per task: CI-mirror (g++-13 Debug+Werror) core 0/0; **build-gcc14** core 0 maybe-uninitialized/0 errors (the GitHub CI compiler — standing gate since the 7c ship failure); maxset core 0 errors. GitHub CI green after ship.
- Identity BOTH modes BITWISE vs `.superpowers/sdd/phase1-identity-reference.txt` — 8a rewrites arithmetic-order-preserving loops; any deviation is a bug, not a legitimate relaxation.
- Physics pattern: `^particle$|cell_system|langevin|brownian|thermostat|npt|dpd|lees_edwards|hybrid|virtual|collision|rigid|bond` (+ engine/rotation for the VV-rotation kernel: `rotation|engine`).
- RNG discipline: thermostat kernels must produce IDENTICAL noise per particle (Philox keyed on id/counter — the kernel must read the SAME key fields in the SAME order; a reordered RNG call = silent identity break that the gate catches — treat identity as the RNG regression test).
- Full-word names.

## Tasks

### Task 1: Hot-loop audit + kernel-form design note
Zero behavior change. For EACH target loop (VV trans 1/2, VV rot, Langevin trans/rot, Brownian trans/rot): record the exact body (file:line), every accessor call and its column, every branch (propagation mask, fixed flags, rotation flag, external-force/virtual guards), the RNG key derivation, and the OMP/Kokkos policy TODAY. Produce the kernel-form design: hoisted pointer set per kernel, branch handling (keep per-row branches — they predict well; NO masking tricks that change arithmetic), where the propagation-mask dispatch stays (outside, per today's hoisted local). Flag any loop whose body calls back into non-columnar code (bond/constraint hooks?) — those keep the view path. Report to `.superpowers/sdd/phase8a-task-1-report.md`. Gates: green baseline (unit, physics pattern, identity). No commit expected.

### Task 2: Column kernels (the flip is per-loop and immediate — no dormant stage; each kernel is small and identity-gated)
- `ParticleStore` gains const-correct raw column accessors for kernel wiring if not already present (`position_data()`, `velocity_data()`, `force_data()`, `mass_data()`, `propagation_data()`, `ext_flag_data()`, `rotation_data()`, gamma per-config, quaternion/omega/torque for rotation — follow the existing `*_view()` host-view getters; add raw-pointer+stride forms).
- Rewrite each target loop per the Task-1 design note: one commit per loop GROUP (VV translation; VV rotation; Langevin; Brownian) so identity bisects trivially. After EACH group: clean rebuild, unit, identity BOTH modes (bitwise), the group's physics tests.
- The `for_each_local_particle` helper stays for non-hot callers.
- Commits `core: velocity-verlet translation as column kernels (phase 8a)` etc.

### Task 3: Measure, iterate, wrap + ship
- Build the phase-7 tip A/B binary once (`/ssd/weeber/eliminate_Particle_builds/ab-phase7`, commit b1831f89f7, standard Release config).
- CONTROLLER-grade measurement (subagent runs foreground sequential, machine idle, check-load before every block): order-balanced A/B ≥5 reps lj1 + lj-omp + p3m1 vs phase-7 tip — target lj1 ≤0.985 (claw back ≥1.5pp); then the full sequential gate 3 reps vs phase-0 — target lj1 ≤1.05 and ALL configs within the ORIGINAL 5% (p3m-omp adjudicated by A/B if the sequential draw is noisy, per the 7c precedent).
- If lj1 misses after honest kernel work: sub-slot instrumentation (the proven std::chrono method), fix what the data indicts within 8a scope, re-measure; if the residual is attributed OUTSIDE 8a scope (ghost PROPRTS streams → the deferred dedup restructure), STOP and report the attribution — do not scope-creep.
- Validation gauntlet (the 7c list + gcc14 gate); spec adjudication (8a results + budget trajectory); push to `me`; GitHub CI green; ledger+memory; hand off to 8b.
