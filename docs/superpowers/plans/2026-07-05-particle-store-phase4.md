# ParticleStore Migration — Phase 4 Implementation Plan (velocity/omega eviction)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evict velocity and angular velocity (omega) from the `Particle` struct into ParticleStore columns; `commit_particle`'s last per-step vector copy dies and kernels read the velocity column directly.

**Architecture:** Exact instantiation of the proven phase-2/3 template: prep audit → columns+carriers → ghost MOMENTUM value path → flip → checkpoint. Velocity is simpler than position (no fold/shift semantics, no derived quantities) but appears in velocity-dependent kernels (DPD, thermalized bonds) via the pack — those reads move to the store column via the existing `pack_index_to_store_row` translation. The phase-3 fetch-cache snapshot store automatically covers Python `v`/`omega_lab` getters once carriers are serialized (no getter rework needed — verify, don't assume). MOMENTUM ghost updates get the same value-branch + columnar-context treatment as POSITION, including membership in the per-step flags when LB/DPD/thermalized bonds are active.

**Spec:** `docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md` (section 3 phase 4, amendment of 2026-07-05). **Templates:** phase-3 flip `git show e98a165212`; ghost value-branch `git show 00ce7e5864`; columnar contexts + hardenings `git show 2d3f84431e e8859513e2`.

## Global Constraints

- All phase-3 global constraints carry over verbatim (worktree-only, `make -j8`, identity gate lj `bdd2022c…`/p3m `141cc4aa…` at every task, clean-core-rebuild after Particle/ParticleStore header ABI changes, test-noise allowlist, controller-only benchmarks, full-word naming, commit trailer).
- Benchmark budget (spec amendment): 1000ppc configs ≤8%, 4000ppc ≤5%, vs `maintainer/benchmarks/baselines/phase0-baseline.csv` (re-recorded 2026-07-05).
- Test gates MUST include `dpd`, `lb`, `thermalized_bond`-related patterns in addition to the standing `lees_edwards|hybrid_decomposition|virtual_sites` set — velocity-dependent physics is this phase's risk surface.
- The T2/T3.5 hot-loop hoists created `auto vel = p.v();` REFERENCE bindings (v() returns `Utils::Vector3d&` today): post-flip these become proxy hoists automatically ONLY where used with `[]`/`=`/`+=`; the audit must fix template-arithmetic and reference-binding uses (`Utils::Vector3d &v = p.v()`, `auto &`).

---

### Task 1: Preparation audit — velocity/omega call sites

**Files:** grep-driven across src/core, src/script_interface, unit tests.

- [ ] Inventory + fix ALL proxy-hostile uses of `v()` and `omega()` per the phase-3 category catalog (.superpowers/sdd/phase3-task-3-report.md): reference bindings → local copy+write-back or proxy-hoist form; template arithmetic (`get_mi_vector`, `sqr`, `operator*` on Vector) → explicit `Utils::Vector3d(...)`; `auto x = p.v()` value-intent → explicit type. Known hot files: integrators (velocity_verlet, npt ×2, brownian, symplectic_euler, stokesian, steepest_descent — many already hoisted as references; convert those hoists to flip-surviving form), thermostats (langevin_inline reads v), dpd.cpp/dpd kernels (velocity-dependent pair forces — CHECK whether they read v via pack (commit_particle velocity) or live particles), lb/particle_coupling.cpp, rattle.cpp (velocity corrections `p.v() += ...`), rotation.cpp (omega, largely hoisted), virtual_sites (relative: omega/velocity transfer; lb_tracers `p.v() = ...`), analysis/statistics (kinetic energy `p.v().norm2()` etc.), observables traits, galilei (kill_particle_motion, momentum), sd_interface, thermalized bond kernel (reads velocities via aosoa pack — bond_forces_kokkos.hpp), EspressoSystemStandAlone/unit tests.
- [ ] Gates: build; unit tests; `ctest -R "integrator|langevin|brownian|npt|dpd|lb|thermalized|rattle|virtual|galilei|analysis"`; identity IDENTICAL. Commit `core: prepare velocity/omega call sites for proxy accessors (phase 4)`.

### Task 2: Store columns + carriers

**Files:** src/core/particle_store/ParticleStore.{hpp,cpp}, src/core/Particle.hpp, unit tests.

- [ ] Columns `m_velocity` (`DualView<double*[3], StateVectorLayout>`) and `m_angular_velocity` (ESPRESSO_ROTATION) + old-generation spares + accessors `velocity_reference/velocity_value/velocity_view` and `angular_velocity_*`; wired into begin_rebuild swap/grow, assign_row preserve_or_seed, finish_rebuild, release_columns.
- [ ] DORMANT carriers on Particle (`m_migration_velocity`, `m_migration_angular_velocity` + `detached_*()`/`migration_*()` getters reading struct fields pre-flip; NOT serialized yet).
- [ ] Unit tests: preservation-by-old-row + carrier seeding + new-row zero defaults for both columns (extend the existing state-column tests).
- [ ] Gates: unit tests; identity. Commit `core: ParticleStore velocity columns (phase 4)`.

### Task 3: Ghost MOMENTUM value path + columnar context

**Files:** src/core/ghosts.cpp (+ ghosts.hpp if the context caches grow).

- [ ] serialize_and_reduce MOMENTUM branch → direction-first value-based (LOAD always writes into the particle; enumerate all four (policy, direction) combinations — the phase-2 lesson), routed through a MomentumRowContext (v+omega raw pointers, hoisted per communication) when store clean & attached, accessor fallback otherwise. Mirror the POSITION branch's structure exactly (a85bf5b105 + 2d3f84431e + e8859513e2 incl. empty-store guard and debug canary).
- [ ] calc_transmit_size: MOMENTUM added to the mask-out scheme (constant 3 or 6 doubles).
- [ ] Columnar whole-communication fast path for `data_parts == GHOSTTRANS_MOMENTUM`? Only if such comms exist (grep callers — rattle does `ghosts_update(POSITION|MOMENTUM)`; per-step LB/DPD flags produce POSITION|PROPERTIES|MOMENTUM mixes). The mixed-parts row-context path is the one that matters; add pack/unpack_momentum_range + LOCL branch ONLY if a pure-MOMENTUM comm exists; otherwise the row-context in serialize_and_reduce suffices — decide from the grep and document.
- [ ] Gates: unit; `ctest -R "dpd|lb|thermalized|rattle|lees_edwards|hybrid|^particle$"` incl. 4-rank; identity. Commit `core: ghost momentum exchange via explicit field values (phase 4)`.

### Task 4: THE FLIP

**Files:** src/core/Particle.hpp (ParticleMomentum sub-struct deleted; serialize rewired with carriers), src/core/short_range_cabana.hpp (commit_particle loses velocity — scalars only remain; pack's velocity view deleted; kernels using velocity read the store column via `row(i)` translation — thermalized bond kernel in bond_forces_kokkos.hpp, any DPD kernel path), src/core/cell_system/CellStructure.{hpp,cpp} (velocity view plumbing to kernels), sweep files, unit-test fixtures.

- [ ] Accessors: `v()` non-const → VectorReference, const → value; `omega()` likewise (ROTATION); assert-attached with carrier fallback for detached (phase-3 adjudicated pattern).
- [ ] Particle::serialize: `ar & m;` leaves; velocity carriers added (SAVE fills from detached getters — column when attached).
- [ ] Verify the fetch-cache snapshot store covers Python `v`/`omega_lab` getters post-flip (they read cached particles via accessors → snapshot columns): run the multi-rank particle tests and a manual 4-rank probe reading p.v of remote particles.
- [ ] Sync-point audit: velocity readers outside existing sync coverage (galilei momentum, analysis momentum/kinetic — most already synced from phase 3; grep for new gaps).
- [ ] Clean core rebuild; full battery `ctest -R "particle|cell_system|integrator|langevin|brownian|npt|dpd|lb|thermalized|rotation|coulomb|lees_edwards|virtual|collision|cluster|analysis|hybrid|rattle|engine|observable|accumulator|galilei|stokesian"`; identity BOTH modes IDENTICAL. Commit flip.

### Task 5: Checkpoint + ship

- [ ] Guard; unit; hybrid until-fail:3; full check_python_skip_long + check_python; identity; HDF5 side-build if h5md-adjacent files changed.
- [ ] CONTROLLER: benchmark gate 3 reps vs baseline — thresholds per amendment (1000ppc ≤8%, 4000ppc ≤5%); commit_particle's death should IMPROVE small-N configs; report the full table and the delta vs the phase-3 numbers (lj1 1.058, lj4 1.057, p3m4 1.048, p3m1 0.981, lj-omp 1.042, p3m-omp 1.010).
- [ ] Spec phase-4 adjudication notes if any deviation; ledger; push to `me`.

## Self-Review Notes
- Velocity-dependent kernels are the phase's unique risk: thermalized bonds read BOTH partners' velocities (pack today), DPD reads relative velocities — the flip must keep kernel index-space discipline (velocity reads by `row(i)`, force writes by pack index `i`).
- Python getter path deliberately reuses phase-3 snapshot-store infrastructure — Task 4 verifies rather than builds.
- MOMENTUM per-step communication only occurs with LB/DPD/thermalized bonds active — the lj/p3m gates don't exercise it per step; correctness comes from the dpd/lb/thermalized test patterns.
