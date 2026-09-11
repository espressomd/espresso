# ParticleStore Migration — Phase 7e Implementation Plan (id→row map; held-pointer sites)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Retire `CellStructure::m_particle_index` (the id→`Particle*` map and its view pool) in favor of the id→row map; `get_local_particle` resolves id→row and returns a view; the latently fragile held-pointer sites (collision detection, bond handler) are rewritten to re-resolve by id/row with debug generation guards — making them safe under the 7c permutation resort that follows.

**Architecture:** Third phase-7 sub-project (order 7a→7b→**7e**→7c→7d; exploration `.superpowers/sdd/phase7-exploration.md` topic 6 is the basis; 7a+7b are SHIPPED at head `1d74f72f0c`). `m_id_to_index` already exists (`Kokkos::View<int*>`, written as pack index == store row on the local prefix, real rows on the ghost tail) — it IS the id→row map; this sub-project makes it the single source of id resolution, kills the view pool + `m_particle_index`, and retires `m_pack_index_to_store_row`/`m_unique_particles` where pack==row makes them identity. 7e lands BEFORE 7c so the held-pointer sites are already row-epoch-safe when rows start physically moving.

**Spec:** section 3 phase 7 (id→index map replaces `m_particle_index` — listed under phase 5 in the table but adjudicated into 7e; spec adjudication note at ship). Perf: expected neutral-to-positive (view-pool refresh per rebuild dies); gate recorded, ship criteria = amended budget (1000ppc ≤12%, 4000ppc ≤5%).

## Global Constraints

- All standing constraints (worktree; `make -j8`; controller benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`).
- MANDATORY before every commit: `maintainer/format/clang-format.sh` (+ autopep8 for touched .py) — CI refuses to run otherwise.
- Stale-ABI rule (11+ incidents): header change → clean-rebuild `espresso_core` + `unit_tests_executables` everywhere tested; NFS `compiler_depend` NUL → `cmake --fresh` regen (Gate-7 lesson: use absolute paths, not the symlink, for build-maxset cmake regen).
- CI gates per header-touching task: CI-mirror core 0/0; maxset core 0 errors. GitHub CI green after ship.
- Identity BOTH modes bitwise vs `.superpowers/sdd/phase1-identity-reference.txt`.
- Physics gates: collision_detection, virtual_sites, bond_breakage suites are THE primary gates (the held-pointer sites); plus the standing pattern `^particle$|cell_system|collision|virtual|lees_edwards|hybrid|dpd|reaction|rigid|bond|exclusion`.
- `get_local_particle` nullptr contract preserved: absent id → nullptr (or a falsy handle) exactly as today; every caller null-checks.
- Full-word names.

## Tasks

### Task 1: id→row map audit + gap check
Zero behavior change. (a) Verify `m_id_to_index` coverage: written for ALL ids (locals + ghosts) at every rebuild? What happens for ids absent on this rank (sentinel)? Ghost dedup: which row does a multiply-imaged ghost id map to (the pack dedup row — record the contract)? (b) Enumerate every `get_local_particle`/`m_particle_index` consumer at current HEAD (exploration topic 6a is stale — 7a/7b changed sites); classify read-and-drop vs held-across-mutation. (c) The view pool: who holds returned `Particle*` beyond one topology epoch? (d) `m_pack_index_to_store_row` + `m_unique_particles`: current consumers; confirm pack index == store row holds everywhere now (the 7a assert), i.e. both are removable. (e) Record the collision-detection queue lifecycle (queue filled during force loop, processed after — where exactly pointers/ids cross a resort boundary). Report to `.superpowers/sdd/phase7e-task-1-report.md`; fix nothing. Gates: green baseline (unit, primary physics pattern, identity). No commit expected.

### Task 2: THE SWITCH — id→row resolution + view-pool death
- `get_local_particle(id)`: resolve via `m_id_to_index` (sentinel → nullptr-equivalent). Return form: keep `Particle*` API by returning a pointer into a SMALL stable id-keyed view cache ONLY if callers require pointer identity across calls — the Task-1 audit decides; PREFER returning by value/optional<Particle> where the caller surface is small enough to convert (views are 16 bytes — by-value is now cheap; convert read-and-drop callers mechanically). Held-across-mutation callers (collision Bind/Glue, bond handler chains) get explicit re-resolution: capture ids, re-resolve after any `add_particle`/topology change, with the 7a `ParticleStoreGuard` generation assert at each re-use point (adopt it — the helper exists).
- Delete: the 7a view pool; `m_particle_index` and its update/rebuild choreography; `m_pack_index_to_store_row` + `m_unique_particles` (pack==row identity — keep the debug assert that proves it at rebuild).
- `m_id_to_index` maintenance moves fully into the store rebuild (single write site; particle_node's map erase interplay per the 7a M3 note — verify remove_particle leaves no stale id entry readable).
- Unit tests: id→row resolution (present local, present ghost, absent id, after add/remove/resort); a collision-detection-shaped regression test: queue a bond-creation across an add_particle (VS creation) and verify the re-resolved partner is correct (the BindAtPointOfCollision pattern — model on existing collision tests; must exercise the re-resolution path).
- Gates: unit ALL; PRIMARY physics (collision/virtual/bond_breakage until-fail:2) + standing pattern; hybrid until-fail:3; identity BOTH modes; CI-mirror; maxset; format scripts. Commit `core: resolve particles by id to store rows; retire the particle pointer index (phase 7e)`.

### Task 3: Wrap + ship
Validation gauntlet (guard, unit, hybrid ×3, full check_python skip-long, identity, CI-mirror full, maxset, HDF5, collision/virtual/bond stress ×2, checkpoint round-trip); CONTROLLER benchmark 3 reps at amended budget; spec adjudication note (id→row map delivered in 7e); push to `me`; GitHub CI green; ledger+memory; hand off the 7c remaining-work list (permutation resort + (offset,count) collapse + the M2/M3 wire cleanups from the 7b review).
