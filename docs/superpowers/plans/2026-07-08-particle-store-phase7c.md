# ParticleStore Migration — Phase 7c Implementation Plan (permutation resort; cells become ranges)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resort becomes a pure column permutation (counting-sort by cell → old_row→new_row permutation → per-column permute kernels), and cells collapse from `Utils::Bag<int>` row bags to `(offset, count)` ranges over the store — the spec's phase-7 storage end state on the CPU.

**Architecture:** Fourth phase-7 sub-project (7a/7b/7e SHIPPED, head `5f1dc4be20`+docs; exploration `.superpowers/sdd/phase7-exploration.md` topic 3 is the basis). The rebuild machinery is ~90% there: `begin_rebuild` double-buffers via generation swap; the per-particle `assign_row` (branchy preserve-or-seed across ~40 columns) is replaced by per-column permute loops (`col_new[new_row] = col_old[perm[new_row]]` — contiguous, vectorizable); ragged sidecars move-on-permute (machinery exists). Ghost rows are NOT permuted (rebuilt fresh by ghost exchange). After the collapse, `Cell::particles()` iterates a contiguous row range (`RowParticleRange` over `[offset, offset+count)` — the CellRows bag dies), mutation windows use the staging store (7b) + pending-removal marks, and the physical order is (re)established only at the permutation rebuild. `m_pack_index_to_store_row`/`m_unique_particles` die IF the pack's ghost dedup can be expressed on rows (Task 1 adjudicates; keep if blocked, with a dated note).

**Spec:** section 3 phase 7 (resort-as-permutation + ranges). Perf: THE banked win — per-column permute replaces assign_row; the phase-7 checkpoint (after 7d) re-evaluates the 1000ppc budget (lj1 currently 1.084; the original 5% is within reach). Gate at ship: amended budget (1000ppc ≤12%, 4000ppc ≤5%) + neutrality-or-better vs the 7e gate CSV (`.superpowers/sdd/phase7e-gate.csv`).

## Global Constraints

- All standing constraints (worktree; `make -j8`; controller benchmarks; commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`).
- MANDATORY before every commit: `maintainer/format/clang-format.sh` on changed files.
- Stale-ABI rule: header change → clean-rebuild `espresso_core` + `unit_tests_executables` everywhere tested; NFS → `cmake --fresh` absolute paths.
- CI gates per header-touching task: CI-mirror core 0/0; maxset core 0 errors. GitHub CI green after ship.
- Identity BOTH modes bitwise vs `.superpowers/sdd/phase1-identity-reference.txt` — the permutation must reproduce the EXACT row order the current rebuild produces (cell-traversal order, surviving-then-staged per cell); this is a pure storage-move phase.
- The 7a Verlet-row generation guard is the safety net for stale rows — extend its coverage if the permutation introduces new consume windows; the physics batteries + collision/virtual/bond_breakage until-fail:2 are the real stale-row gates (identity is blind to them).
- Physics pattern: `^particle$|cell_system|collision|virtual|lees_edwards|hybrid|dpd|reaction|rigid|bond|exclusion` + hybrid until-fail:3 + checkpoint round-trip.
- Full-word names.

## Tasks

### Task 1: Audit + dedup adjudication
Zero behavior change. (a) Record the exact current rebuild order contract (per-cell: surviving rows in bag order, then staged in push order; cells in local-span order, then ghost cells) with file:line — the permutation must reproduce it byte-for-byte. (b) Enumerate every mutation window between rebuilds (insert/extract/erase paths post-7b) and their cell-order effects; design note: how (offset,count) cells express "row removed mid-epoch" (tombstone mark + skip in iteration? immediate compaction is a rebuild — the answer shapes Task 3; today's Bag<int> swap-removes, so within-cell ORDER CHURN on removal is the existing contract — ranges cannot swap-remove, adjudicate: pending-removal list on CellStructure consulted by iteration, resolved at rebuild). (c) Ghost-dedup analysis: can `m_pack_index_to_store_row`/`m_unique_particles` die (pack iterates rows directly; dedup on rows)? Adjudicate with evidence. (d) Inventory RowParticleRange/CellRows consumers for the range collapse. (e) Carries check: 7b-review M2 (ghosts.cpp sizer remnant → RATTLE size arithmetic) + M3 (MigrationPack header ids), 7e-review M1 (stale comment CellStructure.cpp:1067) — fold the trivial ones into Task 2/3 commits. Report to `.superpowers/sdd/phase7c-task-1-report.md`. Gates: green baseline. No commit expected.

### Task 2: Permutation machinery (dormant)
- `ParticleStore::permute_rebuild(std::span<int const> permutation, std::size_t n_local, std::size_t n_ghost)` (or equivalent shape per Task-1 findings): generation swap; per-column permute loops for every column (build from the assign_row column list — three-way sync comment updates: assign_row ↔ copy_row ↔ permute ↔ MigrationPack); ragged sidecars + POD sidecars move-by-permutation; ghost tail allocated fresh (seeded defaults). Unit tests: permutation correctness on a maximally-populated store (identity perm, reversal, random perm — field-for-field vs a reference copy; ragged contents move intact; the maximal-population helpers from 7b Task 2 are reusable).
- Permutation BUILDER on CellStructure (dormant): walk cells in the current rebuild order producing perm[] + the future (offset,count) per cell; ADDITIONAL_CHECKS cross-verify against the live assign_row rebuild result (row-for-row id equality) — run both paths in debug for one release cycle.
- Gates: unit ALL; `^particle$|cell_system`; identity BOTH modes (dormant); CI-mirror; maxset; format. Commit `core: column permutation rebuild machinery (phase 7c)`.

### Task 3: THE FLIP — permutation resort + range collapse
- `ensure_particle_store_synchronized` uses the permutation path: builder → permute_rebuild → per-cell (offset,count) written back; staged insertions land via copy_row into their target range positions (order contract preserved); `Cell` drops the CellRows bag for `(offset, count)`; `RowParticleRange` iterates the contiguous range (iterator simplifies — no row array indirection); mutation windows per the Task-1 adjudication (pending-removal marks consulted by iteration between rebuilds).
- `assign_row` retires from the rebuild hot path (keep only if the staging/fetch stores still need it — check; copy_row may cover them).
- Ghost dedup per Task-1 adjudication: kill `m_pack_index_to_store_row`/`m_unique_particles` if unblocked.
- Fold in: 7b M2/M3, 7e M1 carries (if not already).
- Gates: unit ALL; FULL physics pattern + collision/virtual/bond_breakage until-fail:2 + hybrid until-fail:3 + checkpoint round-trip; identity BOTH modes; CI-mirror FULL; maxset core+unit; guard script (update for the range choke points); format.
- Commit `core: resort as column permutation; cells hold row ranges (phase 7c)`.

### Task 4: Wrap + ship
Validation gauntlet (the 7e list, .superpowers/sdd/phase7e-task-3-validation.md has commands); CONTROLLER benchmark 3 reps: amended budget vs phase-0 AND vs `.superpowers/sdd/phase7e-gate.csv` (expect improvement on 1000ppc; investigate ANY config >1.02 vs 7e with order-balanced A/B before shipping); spec adjudication (range collapse delivered; dedup outcome; budget trajectory note for the phase-7 checkpoint); push to `me`; GitHub CI green; ledger+memory; hand off 7d (GPU shim) + phase-7 checkpoint remaining list.
