# ParticleStore Migration — Phase 7b Implementation Plan (per-field migration pack; carriers die)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace every whole-Particle boost-serialization path (global-resort migration, head-node fetch cache) with per-field packed transfers sourced/sunk directly from ParticleStore columns; delete `Particle::serialize`, ALL `m_migration_*` carriers, the dual-role bond/exclusion members, and `ParticleLocal` — `Particle` becomes a true two-word view `{store*, row}`.

**Architecture:** Second phase-7 sub-project (exploration `.superpowers/sdd/phase7-exploration.md` topic 4 is this plan's basis; 7a landed at commit `4772a9636c`, UNSHIPPED — 7b stacks on it, ship happens as a batch at the 7c gate). The remaining whole-Particle boost users are exactly: (i) `AtomDecomposition::resort` all_to_all, (ii) `RegularDecomposition::exchange_neighbors`, (iii) the head-node fetch cache (`particle_node.cpp` send/recv + gatherv) — ghost exchange is already per-field (phases 3-6). The GHOSTTRANS field grouping partitions the full field set (topic 4 "envelope scope"); migration reuses it with ragged legs for bonds/exclusions. KEY DESIGN ELEMENT (replaces 7a's carrier-based staging): with carriers dead there is no such thing as a detached data-carrying Particle — the 7a staging buffer (`std::vector<Particle>`) and `snapshot_row` are replaced by a **staging ParticleStore** (a second, small store on CellStructure — the `fetch_cache_store` pattern proves the shape): extract = copy row store→staging (per-column, reusing `assign_row` machinery), wire = pack fields from staging rows (or directly from live rows), receive = unpack into staging rows, rebuild = commit staging rows into the main store. `mpiio.cpp`'s detached-carrier vector (T1 audit) moves to the same mechanism.

**Spec:** section 3 phase 7 (per-field pack half). Perf: envelope death is a BANKED gain (exploration topic 10) — sizeof(Particle) → 2 words, no boost archive per migrating particle, no SAVE-time carrier-sync gather. The lj1/lj4 resort residual (whole-store rebuild column copies) is 7c's job — do NOT attack it here. Gate measured at the end but SHIPPING remains deferred to the 7c batch gate (amended budget: 1000ppc ≤12%, 4000ppc ≤5%).

## Global Constraints

- All standing constraints (worktree; `make -j8`; controller benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`).
- MANDATORY before every commit: `maintainer/format/clang-format.sh` (+ autopep8 for .py) on changed files — CI refuses to run otherwise.
- Stale-ABI rule (11 incidents): Particle/ParticleStore/Cell header change → clean-rebuild `espresso_core` + `unit_tests_executables` in EVERY build dir tested; NFS compiler_depend NUL corruption → `cmake .` regen.
- CI gates every header-touching task: CI-mirror core (Debug+Werror) 0/0; maxset core 0 errors.
- Identity BOTH modes bitwise vs `.superpowers/sdd/phase1-identity-reference.txt` (wire-encoding change only — values must be identical).
- Physics gates: the migration canaries are `hybrid_decomposition` (until-fail:3, 4 ranks) and the lees_edwards/virtual_sites migration patterns (P2 lesson); reaction + collision + rigid_bond for the mutation surfaces; full pattern `^particle$|cell_system|collision|virtual|lees_edwards|hybrid|dpd|reaction|rigid|bond|exclusion`.
- The wire format changes; the VALUES must not. Multi-round `exchange_neighbors` (up to grid-sum−3 rounds) must deliver identical final state — the per-field pack must handle a particle migrating through INTERMEDIATE ranks (staged on rank B in round k, re-packed toward rank C in round k+1).
- `reuse_forces`/`run(0)` semantics: forces of migrated particles must survive migration exactly as the carriers preserved them (per-field pack includes the FORCE group).
- Checkpoint format untouched (goes through Python ParticleHandle — verified topic 4).
- Full-word names.

## Tasks

### Task 1: Preparation audit
Zero behavior change. (a) Verify/extend the topic-4 envelope-user inventory at current HEAD (7a moved things): every `boost` archive touchpoint of Particle, every `snapshot_row`/staging-buffer consumer 7a added, `mpiio.cpp`'s vector, the sizer `Particle{}` in ghosts.cpp, `PartCfg`. (b) Field-group completeness table: enumerate `Particle::serialize` legs field-by-field (every ifdef) against the GHOSTTRANS groups + ragged legs — the pack must cover EXACTLY this set; a missed field = silent migration data loss (phase-5 review class). (c) Fetch-cache call graph: `get_particle_data`/`prefetch`/`mpi_get_particles` — which fields callers actually consume (all, via ParticleHandle). (d) Record the multi-round exchange choreography with file:line. Report to `.superpowers/sdd/phase7b-task-1-report.md`; fix nothing unless it blocks Task 2. Gates: green baseline (unit, physics pattern, identity). Commit only if fixes needed.

### Task 2: Staging store + per-field pack machinery (dormant)
- `CellStructure` staging `ParticleStore` (small, lazy-sized) + row-level transfer helpers on ParticleStore: `copy_row(ParticleStore const &source, int source_row, int destination_row)` (every column + sidecars, ifdef-complete — build it from the assign_row column list; unit-test field completeness by round-tripping a maximally-populated particle) — this is the machinery `extract` and `receive` share.
- Per-field pack/unpack free functions (new file, e.g. `src/core/particle_store/MigrationPack.hpp/.cpp`): `pack_rows(store, rows, buffer)` / `unpack_rows(store, first_row, buffer)` covering the FULL field set via the GHOSTTRANS grouping + ragged bond/exclusion legs (length-prefixed, same encoding as the ghost bondbuf) + id-list header. Buffer = `std::vector<char>` with explicit little-endian-native POD memcpy (matches ghost wire practice). Include a versioned size-arithmetic function replacing the sizer-Particle (constant per ifdef config + ragged actuals).
- Unit tests: pack→unpack round-trip identity on a maximally-populated store (every ifdef field, non-empty bonds/exclusions, all POD sidecars); staging-store copy_row completeness; empty/edge cases (zero particles, bond-free).
- DORMANT: production still ships the boost envelope. Gates: unit ALL, `^particle$|cell_system`, identity, CI-mirror, maxset, format scripts. Commit `core: migration staging store and per-field pack machinery (phase 7b)`.

### Task 3: FLIP the migration and fetch paths
- `RegularDecomposition::exchange_neighbors`: ParticleList send buffers → per-field packed buffers (extract = copy_row live→staging + remove from cell + dirty; pack from staging; receive = unpack into staging; `move_if_local` commits staging rows into target cells' pending-commit path; multi-round: staged non-local rows re-pack next round). `AtomDecomposition::resort`: same via all_to_all on packed buffers.
- Fetch cache: owner packs the single row per-field → head unpacks into `fetch_cache_store` row → hands out the view (API unchanged: `get_particle_data` returns `Particle const&`). `mpi_get_particles` gatherv → per-field bulk pack. `PartCfg` copies become fetch-store-backed rows (exploration topic 1 recommendation).
- `mpiio.cpp`: staging-store based read/write.
- 7a's `snapshot_row`/Particle-staging replaced by staging-store rows throughout; ghost `resize_ghost_storage` staging unaffected (ghost rows carry no data until ghost exchange).
- Gates: unit ALL; physics battery FULL pattern; hybrid until-fail:3; particle round-trip (slice/property/checkpoint); identity BOTH modes; CI-mirror; maxset; format. Commit `core: per-field particle migration and fetch (phase 7b)`.

### Task 4: Kill the envelope — Particle becomes two words
- Delete: `Particle::serialize` + boost includes, ALL `m_migration_*` carriers, dual-role bond/exclusion members, `sync_migration_carriers`/`detached_*` choreography, `ParticleLocal l` (`is_ghost()` → `store_row() >= n_local` via the store; `set_ghost` callers — the ghost resize path sets rows, verify the single caller), the sizer `Particle{}` (per-field size arithmetic from Task 2), `snapshot_row`.
- Accessors lose the attached/detached branch (unconditional store access + debug assert attached); constexpr-when-disabled statics stay.
- Unit-test sweep: everything constructing free Particles now attaches to a store first (7a migrated most; finish the rest); serialization tests retire or become pack round-trip tests.
- Expect sizeof(Particle) == pointer + int (+padding); add a static_assert documenting it.
- Gates: unit ALL; physics battery FULL; hybrid until-fail:3; identity BOTH modes; CI-mirror FULL build; maxset core + unit; format; guard script.
- Commit `core: delete the particle migration envelope; Particle is a view (phase 7b)`.

### Task 5: Sub-project wrap (NO ship)
Validation gauntlet (guard, unit, hybrid ×3, full check_python skip-long, identity, CI-mirror full, maxset, HDF5 side-build, bond/collision stress ×2); CONTROLLER benchmark gate 3 reps — RECORD results vs phase-0/phase-6 (expect: envelope death helps lj4; resort residual remains until 7c); NO push (batch ship at 7c gate). Ledger+memory; hand off the 7e/7c remaining-work list.
