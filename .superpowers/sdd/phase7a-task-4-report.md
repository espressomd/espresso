# Phase 7a Task 4 — THE FLIP: cells hand out views (report)

**Worktree:** `/tikhome/weeber/es/.claude/worktrees/eliminate_Particle`
**Branch:** `worktree-eliminate_Particle`
**HEAD at start:** `ccbaf7ecea`
**Date:** 2026-07-07
**Author:** Claude (agent)

---

## Summary

`Utils::Bag<Particle>` is removed from `Cell`. A cell now owns:
- `CellRows m_rows` — `Utils::Bag<int>` of committed ParticleStore row indices;
- `std::vector<Particle> m_staged` — not-yet-committed detached particles;
- `ParticleStore *m_store` — the store its rows index into.

`Cell::particles()` returns a `RowParticleRange` (row-range view over `m_rows` +
store); every particle handed out by the cell system is a non-owning
`{store,row}` view. All cell mutation stays behind `CellParticleStorage`, now
operating on rows + staging. Migration/ghost paths keep working via detached
snapshots + staging + rebuild. The transitional `m_particle_index` (id→`Particle*`)
now points into a stable per-CellStructure view pool refreshed on every rebuild.

---

## Staging design (insert / rebuild)

- `CellParticleStorage::insert_particle(Cell&, Particle&&)` appends the incoming
  detached, carrier-laden particle to `cell.staged()` and returns a reference to
  it. The caller marks the store dirty. The particle becomes a committed row
  (columns seeded from its carriers via `assign_row`) only at the next store
  rebuild.
- `ensure_particle_store_synchronized()` (CellStructure.cpp) now, per cell in
  cell-traversal order (local cells then ghost cells):
  1. captures the cell's surviving committed row indices;
  2. clears the row bag;
  3. for each surviving row, hands `assign_row` a `Particle` attached to the
     OLD row so the column data is PRESERVED-by-old-row from the spare columns;
  4. for each staged particle, hands `assign_row` the detached staged particle
     so its row is SEEDED from carriers;
  5. clears `cell.staged()`.
  The store's `n_local`/`n_ghost` are counted as `rows + staged` per cell before
  `begin_rebuild`. Row-assignment order (local prefix `[0,n_local)`, ghost suffix,
  cell order, bag order within a cell) is byte-identical to pre-flip.
- `add_particle` / `add_local_particle` stage + mark dirty + **commit
  immediately** (`ensure_particle_store_synchronized`) and return a stable
  view-pool pointer, preserving the pre-flip contract that a newly added
  particle is live (indexed, iterable, readable/writable) right after the call.
  `set_particle_decomposition` and the migration re-insert paths batch (stage
  all, one rebuild) to avoid O(n^2).

## snapshot_row semantics

`Particle ParticleStore::snapshot_row(int row)`:
- attaches a fresh `Particle` to `(this, row)`;
- calls `Particle::detach_from_store()`, which (a) captures the STRUCTURAL ghost
  flag (`row >= n_local`) into the fallback `l.ghost` member, (b) runs
  `sync_migration_carriers()` — the SAME field set/order as the serialize SAVE
  block, refactored into a private method reused by both — copying every
  column/sidecar value into the migration carriers, then (c) nulls the store
  pointer / sets row -1.
- The result is a self-contained detached value the migration/resort paths move
  into send buffers or a cell's staging area; it does NOT alias the store.
`extract_row(Cell&, index)` = `snapshot_row(cell.rows()[index])` + swap-with-back
remove of that row index (mirrors the pre-flip Bag erase; order not preserved).

## View-pool contract (`m_particle_index` transitional form)

- `m_view_pool` is a `std::deque<Particle>` of attached views, one per indexed
  particle (locals in row order, then the ghost rows whose id is not owned by a
  local — "locals win"). `std::deque` so refilling never relocates
  already-handed-out element addresses within a generation.
- `m_particle_index` (`std::vector<Particle*>`, id-indexed) points into the pool.
  `get_local_particle(id)` returns the pool pointer — valid between rebuilds,
  exactly the pre-flip pointer-stability contract.
- Rebuilt by `rebuild_particle_index()` (LOCALS only; their id columns are valid
  after the rebuild) inside `ensure_particle_store_synchronized`, then
  `index_ghost_particles()` (ghost tail) AFTER `ghosts_update` has filled the
  ghost id columns (a fresh ghost row carries a default id until then).
- `resort_particles` clears the index, resorts, then commits + reindexes so the
  index is consistent on return (a direct `resort_particles` caller — e.g.
  `check_particle_index` — sees a populated index).
- `set_index_map`'s `m_unique_particles` (`std::vector<Particle*>` consumed by
  the Cabana kernels) is rebuilt to point into the pool (its local prefix ==
  pack order == store rows `[0,n_local)`; ghost tail == deduped ghosts), since it
  can no longer hold pointers into transient cell iteration.

## Migration / ghost choreography

- `RegularDecomposition::resort` / `AtomDecomposition::resort` /
  `HybridDecomposition::resort`: iterate a cell's committed rows BY POSITION;
  `fold_and_reset` writes through a view into the column; `extract_row`
  snapshots + removes the row (swap-with-back → re-examine the swapped-in
  position); the snapshot goes into a `ParticleList`/`std::vector<Particle>`
  send buffer (exempt Bags, unchanged) or another cell's staging via
  `insert_particle`. Received particles are re-inserted via `insert_particle`
  (staging). `CellStructure::resort_particles` wires every cell's store before
  the decomposition runs and commits + reindexes after.
- Ghost communication (`ghosts.cpp`): `GhostCommunication::part_lists` changed
  from `std::vector<ParticleList*>` to `std::vector<Cell*>`. Loops iterate
  `cell->particles()` (views); PARTNUM resize routes through
  `resize_ghost_storage(Cell&, count)` (stages `count` default ghosts).
  `contiguous_store_rows` reads the row range directly from `cell.rows()`.
- **Ordering fix (canary):** `update_ghosts_and_resort_particle` now runs
  `resort → COMMIT (sync) → ghosts_count(PARTNUM) → COMMIT (sync) →
  ghosts_update(DATA) → index_ghosts`. The commit before PARTNUM is required so
  source (local) cells report their full committed row count.
- **Multi-layer ghost forwarding fix:** ghost resize STAGES (deferred commit),
  so within one PARTNUM communicator pass a cell resized as a ghost destination
  must report its pending size when read as a source for a downstream ghost
  layer. `Cell::size()` = `rows + staged`, and all PARTNUM count paths
  (`cell_cell_transfer`, `prepare_send_buffer`) use it.
- **Hybrid internal ghost comms:** `HybridDecomposition::resort` does its own
  PARTNUM/DATA ghost comms after the child resorts (which stage). A
  `set_commit_store` callback (wired by `set_hybrid_decomposition` to
  `ensure_particle_store_synchronized`) commits before/between those comms.

## Cached-view multipass hazard (the systemic pitfall)

`RowParticleIterator::operator*` returns a reference to a view CACHED INSIDE the
iterator, overwritten on `++`. Any code that STORES `&p` from a
`for (auto &p : range)` loop (or an STL algorithm that holds references across
increments) dangles. Fixed at every such site:
- `fetch_particles` (observables): rebuilt from `get_local_particle(id)` stable
  pool pointers.
- `lb_couple_particles`: snapshot coupled views into an owning `std::vector<Particle>`
  (view copies alias the store, so force write-back lands in columns), pass
  pointers into that.
- `dipolar_direct_sum::gather_particle_data`: returns `std::vector<Particle>`
  (view copies) instead of `std::vector<Particle*>`; force/torque/dip_fld
  write-back through them lands in columns.
- `cells.cpp get_interacting_neighbors`: collects neighbor IDS instead of
  Particle pointers.
The verlet consume branch and the parallel `for_each`/`reduce` single-cell paths
build one view per index (thread-safe) or reuse one cached-view iterator per cell.

## Guard-script changes

`maintainer/CI/check_cell_storage_mutations.sh` now flags mutating calls on
`rows()` / `staged()` (and, defensively, the legacy `particles()`):
`insert|erase|clear|resize|push_back|emplace_back|pop_back|emplace`. Exempts
`ParticleListOperations.hpp` (the choke points) and `CellStructure.cpp` (the
store rebuild legitimately clears + refills each cell's row bag). Verified it
still CATCHES an out-of-band `cell.rows().insert(...)` (synthetic tripwire test)
and passes clean on the tree.

## Generation-invariant re-verification (T3 table, re-checked at the flip)

The verlet list holds store rows stamped with `(store, generation)`; a consume
is safe iff no generation bump occurred since the stamp without a verlet rebuild.
Every mark-dirty site and its verlet pairing after the flip:

| Site | Marks dirty | Bumps generation? | Verlet rebuild before next consume? |
|---|---|---|---|
| ctor | yes | on first sync | initial `m_rebuild_verlet_list=true`; safe |
| `remove_particle` | yes | at next resort's sync | removal → resort → flag set |
| `add_local_particle` | yes | IMMEDIATE sync | always followed by a resort before a force calc (controlled placement); debug guard covers |
| `add_particle` | yes | IMMEDIATE sync | sets a resort level → resort → flag set |
| `remove_all_particles` | yes | next sync | followed by resort |
| `ghosts_count` | yes | next sync (inside update_ghosts_and_resort) | only reached after `resort_particles` set the flag |
| `resort_particles` | yes | its own end-of-function sync | sets `m_rebuild_verlet_list=true` in the same function |
| `set_particle_decomposition` | yes | its own sync | full rebuild; next resort sets the flag |

New at the flip: the add paths and `resort_particles`/`set_particle_decomposition`
now bump the generation via an IMMEDIATE sync (not lazily). The verlet consume
branch is only reached from the CPU verlet/collision path, always after a resort
(which sets the flag → BUILD branch restamps). No production path consumes the
verlet list at a generation newer than its stamp. The debug
`ParticleStoreGuard::assert_generation` in the consume branch + the new
`execute_bond_handler` entry guard are the defense-in-depth canaries.

## ImmersedBoundaries dangling-ref hazard

`resolve_bond_partners` returns `Particle*` from `get_local_particles` →
`get_local_particle`, which now returns STABLE view-pool pointers. The
`Particle &p2 = *partners[0]` pattern therefore stays valid for the handler call
duration (no rebuild during a bond kernel). A debug generation guard is added at
`execute_bond_handler` entry, asserting the store generation is unchanged across
each handler call.

## What remains for 7b / 7c / 7e

- **7b:** the migration/ghost WIRE format is still whole-`Particle` boost
  envelopes + the migration carriers; the fetch cache / PartCfg / mpiio still use
  detached-Particle value snapshots. `snapshot_row` and the carrier sync are the
  seam 7b replaces with per-field column packing. `mpiio` add loop and the
  dipolar/lb view-copy buffers are O(n)-per-add / whole-view copies to revisit.
- **7c:** `CellRows` (`Bag<int>`) collapses to a contiguous `(offset,count)`
  range once resort becomes a column permutation; staging + the per-cell row bag
  go away.
- **7e:** `m_particle_index` (`Particle*` + view pool) retires for a pure id→row
  map; `m_unique_particles` dies with the pack; the collision/bond held-pointer
  sites adopt the row/generation model.

---

## Debugging lessons (bugs found + fixed during the flip)

1. **Ghost-index drop on bare sync.** `rebuild_particle_index` initially indexed
   only locals; a bare sync (add_particle commit, forces.cpp) dropped the ghost
   index entries a prior ghosts_update established -> `get_local_particle(ghost_id)`
   returned nullptr -> LB coupling / collision null-deref. Fixed: rebuild also
   indexes already-valid ghost rows (id >= 0); fresh ghosts (id -1) are added by
   `index_ghost_particles` after ghosts_update.
2. **Stale ghosts re-indexed after remove-all.** `remove_all_particles` cleared
   only local cells; stale ghost rows got re-indexed -> a fresh add of the same
   id saw "already exists". Fixed: clear ghost cells too.
3. **Cached-view multipass hazard.** Storing `&p` from a `for (auto &p : range)`
   loop over the row-view range dangles (the cached view is overwritten on ++).
   Hit in fetch_particles (observables), lb_couple_particles, dipolar
   gather_particle_data, cells.get_interacting_neighbors. Fixed via stable
   pool pointers or owning view-copy buffers (see the multipass section above).
3. **PARTNUM ordering / pending-size.** Ghost count must run after locals are
   committed; a cell resized as a ghost dst reports rows+staged (Cell::size())
   for downstream ghost layers in the same PARTNUM pass.
4. **resize_ghost_storage does not mark the store dirty.** The hybrid's internal
   commit-callback therefore has to `mark_particle_store_dirty()` before syncing,
   or the staged ghosts never commit (the bare sync early-returns on a clean
   store) -> MPI_ERR_TRUNCATE in the hybrid ghost DATA transfer. This was the
   hybrid_decomposition (canary) failure.
6. **Staged inserts invisible within the same resort (hybrid ordering).** In
   HybridDecomposition::resort the type-based moves stage particles into their
   target child cells; pre-flip those inserts were immediate and VISIBLE to the
   child resorts that follow. Post-flip they are staged (not committed), so the
   child resorts iterated a different set -> different final placement/order ->
   1-ULP energy drift in `test_resort`. Fixed by committing (commit-callback)
   between the type-moves and the child resorts, restoring pre-flip visibility.

5. **add_particle immediate commit vs teardown.** add_particle/add_local_particle
   commit immediately (so the returned pointer + index are live). A unit test
   that never resets the global System then releases store Kokkos columns after
   Kokkos::finalize -> abort; fixed by resetting the System in that test's main.

## Gates

| Gate | Result |
|---|---|
| C++ unit tests ALL (`build`, `ctest -L unit_test`) | **147/147 passed** |
| Identity `--mode lj` | `bdd2022cb08c8c1c74fe0860ccb53dd0d9833966e97fb2f8f2d8f7232c5d7a69` -- matches reference (BITWISE) |
| Identity `--mode p3m` | `141cc4aa9343884dd57756cd612343cc99016c59bc87b6f575e35b7fa02ee6ce` -- matches reference (BITWISE) |
| Physics battery (full pattern) | **72/74**; the 2 failures are the pre-existing stale-`*_processed.py` build-dir artifacts (sample_reaction_methods_*), which PASS after `find build/testsuite -name '*_processed.py' -delete` -- documented environmental in T1/T2/T3, not code |
| hybrid_decomposition (canary) `--repeat until-fail:3` | **3/3 PASS** (1-core + 4-rank) |
| Particle round-trip (particle_slice\|particle_property\|checkpoint) | **26/26 passed** |
| Guard script | PASS (verified it CATCHES a synthetic `rows().insert()`; clean on tree) |
| CI-mirror core (Debug + -Werror) | **0 errors / 0 warnings** |
| Maxset core (ADDITIONAL_CHECKS) + key maxset unit tests | **0 errors**; 14/14 maxset unit tests pass under runtime ADDITIONAL_CHECKS |

Commit: `b169cd56aa` (`core: cells hold store rows and hand out particle views (phase 7a)`).
