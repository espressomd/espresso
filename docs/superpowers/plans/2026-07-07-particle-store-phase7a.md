# ParticleStore Migration — Phase 7a Implementation Plan (cells hold rows; Particle becomes a view)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cells stop owning Particle objects: `Utils::Bag<Particle>` in cells is replaced by a cell-local container of store ROW INDICES, and every particle handed out by the cell system is a non-owning view (`{store*, row}`) constructed on demand; all cell mutation stays behind the `CellParticleStorage` choke points, now operating on rows.

**Architecture:** First of five phase-7 sub-projects (7a → 7b → 7e → 7c → 7d, exploration report `.superpowers/sdd/phase7-exploration.md` — this plan's basis, read it before each task). ADJUDICATED DEVIATION from the exploration's 7a sketch, to keep 7a shippable without the permutation resort (7c): cells hold `Utils::Bag<int>` row indices (cheap swap-with-back mutation of ints), NOT yet `(offset,count)` ranges — the collapse to contiguous ranges happens in 7c when resort becomes a column permutation. Likewise the accessor attached/detached branch and ALL migration carriers/dual-role members SURVIVE 7a (migration and the fetch cache still ship whole-Particle boost envelopes until 7b); what 7a removes is cell OWNERSHIP. Store rebuild (`ensure_particle_store_synchronized`) keeps working: it traverses cells in order reading each row-index bag, and `assign_row` keeps seeding from the old row via the view's `store_row()`.

**Spec:** section 3 phase 7 (this sub-project delivers the "cells hold store references" half; ranges land in 7c). Perf budget: no perf claim banked at 7a (exploration topic 10); the benchmark gate still runs at the amended budget (1000ppc ≤12%, 4000ppc ≤5%) to catch accidents.

## Global Constraints

- All standing constraints: this worktree only; `make -j8`; controller-only benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Stale-ABI rule (TEN incidents): after ANY change to Particle.hpp/ParticleStore.hpp/Cell.hpp/CellStructure.hpp/pack headers, clean-rebuild `espresso_core` AND `make -C build -j8 unit_tests_executables` in EVERY build dir you test; NFS `compiler_depend.make` NUL corruption → `cmake .` regenerate.
- CI gates at every header-touching task: CI-mirror core compile (build-ci-mirror symlink, Debug+Werror) 0 errors/0 warnings; maxset core compile (build-maxset) 0 errors. GitHub Actions green after the sub-project push.
- Identity gate BOTH modes vs `.superpowers/sdd/phase1-identity-reference.txt` — 7a is a pure storage-ownership move: BITWISE identity must hold.
- Physics gates: full ctest on 1 AND 4 ranks at the ship gate; every task runs at least `^particle$|cell_system|collision|virtual|lees_edwards|hybrid|dpd` (ghost/storage change pattern).
- The transitional invariant `store.id(p.store_row()) == p.id()` (debug builds) must keep holding after every mutation/resort.
- Full-word names over abbreviations.
- CRITICAL KNOWN RISK (exploration topic 6c): `m_verlet_list` holds `std::pair<Particle*,Particle*>` across steps. In 7a, cells no longer own Particles, so THERE ARE NO STABLE Particle* TO POINT AT. Task 3 addresses this surface explicitly BEFORE the flip task; the identity gate is blind to stale-row bugs — the physics battery and a targeted regression test are the real gate.

## Tasks

### Task 1: Preparation audit
Grep-driven, zero behavior change. (a) Enumerate every `Particle*`/`Particle&` that outlives the statement that produced it (m_verlet_list pairs, Cell::m_verlet_list, m_particle_index, m_unique_particles, collision-detection held pointers, virtual-sites reference particles, bond-handler static_vector<Particle*> — exploration topics 6a/6c have the inventory; verify and extend it). Classify: dies-in-7a (must be re-expressed as rows/ids in Task 3) vs survives-7a (rebuilt every use). (b) Enumerate every by-value `Particle` copy/move (Bag internals, extract_particle, MPI buffers, PartCfg, fetch cache — topic 1 table) and mark which 7a touches (cell Bags only) vs 7b (migration/fetch). (c) Find every direct `.particles()` Bag-API use that is NOT iteration/size (mutation outside CellParticleStorage would be a phase-1 guard violation — the guard script should prove none exist). Fix nothing unless it would break under Task 4; write the audit to `.superpowers/sdd/phase7a-task-1-report.md`. Gates: green-baseline run (unit, physics pattern above, identity). Commit only if fixes were needed: `core: prepare particle pointer lifetimes for view particles (phase 7a)`.

### Task 2: Row-bag storage + view-Particle construction machinery (dormant)
- `Cell`/`ParticleList` side: introduce the row-index container type (`Utils::Bag<int>` behind a name like `CellRows`) and a parallel member on `Cell` alongside the existing `ParticleList` (dormant — filled during store rebuild, verified against the Bag contents under ADDITIONAL_CHECKS, not yet read by anyone).
- `ParticleStore` side: `make_view(int row) -> Particle` factory (constructs an attached view Particle bound to this store+row — constructor exists since phase 2 via the attach handle; add the explicit factory + assert row bounds).
- `ParticleIterator`/`ParticleRange` groundwork: a row-range iterator adaptor that yields `Particle` views from a `CellRows` + store (dormant, unit-tested standalone).
- `ensure_particle_store_synchronized` fills the row bags as it assigns rows (one extra int write per particle per rebuild).
- Unit tests: row-bag fill matches Bag traversal order; view factory identity (`store.make_view(r).id() == store.id(r)`); iterator adaptor yields the same sequence as Bag iteration on a built store.
- ABI → clean rebuilds; CI-mirror + maxset. Commit `core: cell row-index storage and particle view machinery (phase 7a)`.

### Task 3: Pointer-lifetime hardening (pre-flip)
Re-express the dies-in-7a pointer surfaces from the Task-1 audit as row/id-based BEFORE cells flip, while both representations exist:
- `m_verlet_list` (CellStructure.hpp:188 and Cell::m_verlet_list): pairs of `Particle*` → pairs of store ROWS, rebuilt on Verlet rebuild as today PLUS invalidated on store-generation change (assert in debug: the list's recorded generation matches the store's before use — the generation counter exists). The pair force kernels resolve rows → views at loop entry (hoisted store pointer, view per row — same cost shape as today's pointer deref).
- `m_unique_particles` (`std::vector<Particle*>`): → rows or eliminate early if the pack path allows (exploration says it dies with the pack in 7e — only re-express if it survives Task 4; do the minimal form).
- REGRESSION TEST (proven protocol): a test that resorts/permutes the store between Verlet build and force use and would read the WRONG particle through a stale entry — must FAIL if the generation invalidation is removed (demonstrate by temporarily disabling it), PASS with it.
- collision-detection/bond-handler held pointers: NOT rewritten here (they re-resolve within one topology epoch and remain valid under 7a because rows don't move between rebuilds; they are 7e's job) — but add the debug generation-assert helper they can adopt.
- Gates: unit + physics battery (lees_edwards/hybrid/virtual/collision emphasized) + identity. Commit `core: verlet list holds store rows instead of particle pointers (phase 7a)`.

### Task 4: THE FLIP — cells hand out views
- `Cell::particles()` returns the row-range view (iterator adaptor over `CellRows` + store) instead of `Bag<Particle>&`; the `ParticleList Bag` leaves `Cell`. `CellParticleStorage` choke points rewritten to row ops: `insert_particle` = append a detached Particle's data to the store?? NO — insertion of NEW particles still goes through the dirty-store rebuild path exactly as today: the incoming (detached, carrier-laden) Particle is stashed in a staging area (a small owning buffer on CellStructure, e.g. `std::vector<Particle>` of pending insertions), the store marked dirty, and the rebuild assigns rows and seeds columns from carriers (assign_row's existing seed branch). `extract_particle` = snapshot the row into a detached carrier-laden Particle (add explicit `Particle ParticleStore::snapshot_row(int row)` using the existing serialize-sync choreography), remove the row index from the cell bag, mark dirty. `erase_particle`/`clear_particles`/`resize_ghost_storage` analogous (ghost resize: stage ghost slots, rebuild assigns ghost rows).
- `ensure_particle_store_synchronized` iterates: cell row bags for surviving particles (preserve via old row) + staging buffers (seed from carriers) in cell-traversal order; staging cleared after rebuild.
- `m_particle_index` (`Particle*`) cannot point into Bags anymore: transitional form = rebuild it to point into... NOTHING stable exists. Interim: replace `get_local_particle` internals with id→row lookup returning a pointer to a per-CellStructure view-pool entry (a stable `std::deque<Particle>` of views refreshed on rebuild — the pointer stays valid between rebuilds, exactly today's stability contract) — this keeps the `Particle*` API for 7e to retire.
- Migration/ghost paths keep working via snapshot/staging (verify RegularDecomposition::exchange_neighbors + AtomDecomposition round-trip: extract→snapshot→wire→staging→rebuild).
- Delete: Bag<Particle> from cells, ParticleListOperations Bag paths, the phase-1 primitives' Particle-moving forms (update the guard script's expectations!).
- Unit tests: all free-Particle-in-cell tests migrate to the store-attach pattern; link_cell/ParticleIterator tests against row ranges.
- Gates: unit ALL; physics battery full pattern + reaction + collision + rigid_bond; particle round-trip (particle_slice|particle_property|checkpoint); identity BOTH modes; CI-mirror + maxset; guard script updated+passing.
- Commit `core: cells hold store rows and hand out particle views (phase 7a)`.

### Task 5: Checkpoint + ship
Guard; unit; hybrid until-fail:3; CI-mirror FULL build; full check_python(+skip_long) on 1 AND 4 ranks; identity; HDF5 side-build; CONTROLLER benchmark gate 3 reps at amended budget (1000ppc ≤12%, 4000ppc ≤5%) with order-balanced A/B vs the phase-6 tip binary (build it to /ssd/weeber/eliminate_Particle_builds/ab-phase6 once) if the sequential gate is within 2pp of a boundary; spec adjudication (document the Bag<int> transitional deviation); push to `me`; verify GitHub CI green; ledger+memory.
