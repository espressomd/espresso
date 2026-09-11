# Array-Based Particle Storage in the ESPResSo Core — Design

Date: 2026-07-03
Status: approved (all sections reviewed during brainstorming)

## Goals

- Replace the `Particle` struct (AoS) as the core's particle storage with
  array-based (SoA) storage that is both SIMD- and GPU-friendly, using Kokkos.
- Keep `Particle` objects in the Python interface; the user-facing API does
  not change.
- Clear separation between **state** (position, quaternion, velocity, mass,
  ...) and **observables** (forces, torques, ...): they live in separate
  containers.

## Decisions made during brainstorming

| Topic | Decision |
|---|---|
| Migration strategy | Incremental via proxy: SoA becomes source of truth; a `Particle`-compatible view/proxy keeps all subsystems compiling; hot paths migrate to direct array access subsystem-by-subsystem. |
| Vector layout | Component-major: one `Kokkos::View<double*[3], LayoutLeft>` per quantity (memory: all x, all y, all z — "x/y/z as separate arrays" with a single handle). |
| Residency | Per-field dual residency (DualView-style host+device pair with per-field modified/sync tracking). `Kokkos::ScatterView` for force/torque accumulation. |
| Ragged data (bonds, exclusions) | Mutation-friendly host sidecar containers as source of truth; flattened CSR-style Kokkos views rebuilt for kernels (existing `LocalBondState` pattern). |
| Store shape | Cell-sorted flat store: locals contiguous sorted by cell, ghosts appended; `Cell` = (offset, count) range; resort = stable permutation of all columns. |
| Dependencies | Storage layer is pure Kokkos; Cabana kept only for neighbor lists (`Cabana::VerletList`, `neighbor_parallel_for`). |
| Success bar per checkpoint | Test suite always green. Perf: at most 5% cumulative regression relative to the phase-0 baseline on LJ and P3M benchmarks at: (1) `--particles_per_core 1000` on 1 and 4 MPI ranks; (2) `--particles_per_core 4000` with 4 OMP threads. |

**Amendment (2026-07-05, user-approved):** for the phase-3 and phase-4
checkpoints only, the `--particles_per_core 1000` configurations (all rank
counts) are allowed up to 8% cumulative regression under the sequential gate
protocol (the `--particles_per_core 4000` configurations stay at 5%).
Same-state interleaved A/B measurements place lj-1rank at ~+3%; the gate's
sequential protocol adds cross-run machine drift on this P/E-core host.
Rationale: measured, bounded residual cost of columnar storage on the
ghost-heavy small-per-rank regime (multi-stream ghost packing, store rebuild
machinery) after seven verified optimizations; it amortizes with per-rank
particle count and is scheduled to shrink structurally in phase 4
(commit_particle removal) and phase 5 (columnar PROPERTIES ghost exchange).
The 5% budget for all configurations is re-tightened at the phase-5
checkpoint. The phase-0 baseline was re-recorded 2026-07-05 from a rebuilt
phase-0 binary (machine-state drift ~3-4%; see baselines/machine-info).

**Amendment (2026-07-06, user-approved):** the phase-5 checkpoint re-tightening
did not hold for the `--particles_per_core 1000` LJ configurations. Final gate
(3 reps, quiet machine): lj-1rank 1.113, lj-4rank 1.103, lj-omp 1.044,
p3m-1rank 0.991, p3m-4rank 1.051, p3m-omp 1.027. A perf-recovery round
(pack-contiguous type cache, solver-gated packed charge/dipm, integrator
accessor hoists) recovered p3m-4rank fully (interleaved A/B: 4.4% FASTER than
phase 4; its 1.051 gate reading is cross-run baseline drift) but left an
lj residual that sub-slot wall-timer instrumentation shows is diffuse SoA
multi-stream cost (per-column PROPRTS ghost serialization, integrator
half-step, pair/pack cache streams) with no single recoverable hot spot —
re-confirming the phase-3.5 finding of an inherent ~6-11% columnar penalty at
1000 particles/rank on this host. Amended bar from phase 5 onward: the
`--particles_per_core 1000` configurations are allowed up to 12% cumulative
regression; the `--particles_per_core 4000` configurations stay at 5%. The
1000-ppc budget is re-evaluated at the phase-7 checkpoint, where the Particle
struct, migration carriers, and accessor attached/detached branches are
removed.

**Phase-6 checkpoint (2026-07-07, PASSED amended budget):** lj-1rank 1.112,
lj-4rank 1.092, lj-omp 1.048, p3m-1rank 0.996, p3m-4rank 1.036, p3m-omp 1.028.
Phase 6 verified perf-neutral vs phase 5 by order-balanced interleaved A/B
(lj1 1.000, lj-omp 1.003, p3m1 1.000) after one recovery fix: the has-exclusion
pack flag became rebuild-cadence and the per-step pack-commit loop — a pure
no-op once phases 5/6 store-aliased every per-step field — was deleted.
Metrology note: sub-1% A/B on this host requires ORDER-BALANCED sampling
(alternate which binary runs first per repetition; the second runner carries a
consistent thermal bias).

**Phase-7a/7b checkpoint (2026-07-07, PASSED amended budget):** lj-1rank 1.082,
lj-4rank 1.014, lj-omp 0.996, p3m-1rank 1.033, p3m-4rank 1.016, p3m-omp 1.014
vs phase-0 — every configuration also inside the ORIGINAL 5% budget except
lj-1rank (covered by the 12% amendment). Adjudicated transitional deviation
(7a): cells hold `Utils::Bag<int>` row indices instead of the spec's
`(offset,count)` ranges until the phase-7c permutation resort lands (ranges
cannot support cheap mutation without it). Phase-7b deleted the whole-Particle
boost migration envelope and every migration carrier: migration and the
head-node fetch path ship per-field packed buffers (MigrationPack), and
`Particle` is now a 16-byte non-owning view `{ParticleStore*, row}`. The
envelope death recovered the phase-7a resort residual and more (lj-4rank
0.824 vs 7a). Remaining 7 sub-projects: 7e (id→row map), 7c (permutation
resort + range collapse), 7d (GPU compile-safe shim).

**Phase-7e checkpoint (2026-07-08, PASSED):** lj-1rank 1.084, lj-4rank 1.020,
lj-omp 0.998, p3m-1rank 1.038, p3m-4rank 1.004, p3m-omp 1.029 vs phase-0;
neutral vs 7b (≤1.015 all configs). The id→`Particle*` index and the 7a view
pool are retired: `get_local_particle` resolves through a store-rebuild-cadence
id→row map and returns views by value; the collision handlers re-resolve by id
after topology changes (a latent GlueToSurface stale-reference bug was found
and fixed with a deterministic regression test). Adjudication: the spec's
"id→index map replaces m_particle_index" (listed at phase 5) was delivered
here; `m_pack_index_to_store_row`/`m_unique_particles` survive until 7c (row
identity holds only on the local prefix).

**Phase-7c checkpoint (2026-07-08, PASSED):** resort is a column permutation
(counting-sort by cell → per-column permute) and cells hold `(offset,count)`
ranges — the spec's phase-7 CPU storage end state; the 7a `Bag<int>`
transitional deviation is retired. lj-1rank 1.066, lj-4rank 1.002, lj-omp
0.990, p3m-1rank 1.036, p3m-4rank 1.010 vs phase-0; p3m-omp read 1.086 on the
sequential gate but order-balanced A/B (6 reps) against the 7e binary gives
0.995 with ±8% intra-binary P3M auto-tuning spread — adjudicated as tuning
variance, config passes (A/B min 1.003 vs phase-0). Mid-epoch removals are
pending-removal marks resolved at rebuild; the removal-order contract change
is invisible to removal-free histories (identity stayed bitwise). One
regression was caught at the gauntlet and fixed: hybrid decomposition's nested
type-migration loops re-staged already-migrated particles once per regular
cell (125×) — fixed with pending-removal guards. The pack ghost-tail
translation maps survive (deduped tail is non-contiguous); they die when ghost
dedup is restructured. Remaining: 7d (GPU compile-safe shim), phase-7
checkpoint budget re-evaluation.

**Phase-8a checkpoint (2026-07-08, PASSED amended budget):** all 8
integrator/thermostat hot loops (velocity-Verlet translation+rotation,
Langevin, Brownian) are direct column kernels via a shared row launcher;
identity stayed BITWISE per kernel group. Perf verdict (order-balanced A/B +
sequential gate): 8a is ≈neutral vs the phase-7 tip (the arithmetic was
already fast; a 2-D-view-subscript regression found by slot instrumentation
was fixed by per-row VectorReference bundles); vs phase-0: lj-1rank 1.086,
lj-4rank 1.014, lj-omp 1.007, p3m-1rank 0.985, p3m-4rank 0.971, p3m-omp
1.033 — five of six configs inside the ORIGINAL 5% budget. The lj-1rank
residual is attributed (instrumented) to the phase-7 whole-store rebuild per
resort, not to accessor proxies — the 12% amendment for 1000ppc stays until
the phase-8 final checkpoint. Doxygen reference cleanup (51→0 warnings)
joined the standing gates alongside gcc-14.

**Phase-8 FINAL checkpoint (2026-07-08, PASSED THE ORIGINAL 5% BUDGET — the
2026-07-05/06 amendments are RETIRED):** lj-1rank 1.049, lj-4rank 0.976,
lj-omp 0.982, p3m-1rank 0.981, p3m-4rank 0.947, p3m-omp 0.993 vs the phase-0
AoS baseline — four of six configurations are FASTER than the original
struct-of-arrays code. The last ~4pp on lj-1rank came from hardware-counter
profiling (perf, enabled on the host late in phase 8): (i) gcc had stopped
inlining the Langevin friction call tree inside init_forces_and_thermostat
(per-particle PLT calls, 15% of core time — fixed with [[gnu::flatten]]);
(ii) the phase-7c pending-removal mask was consulted per element on every
iteration path (fixed with an O(1) has_pending_removals fast path); (iii) the
Verlet-build kernels paid per-candidate Particle-accessor derefs for id and
position (fixed with raw column pointers hoisted per work item). All fixes
bitwise-identity-preserving. Upstream finding (NOT this branch's): p3m_gpu
tuning with zero particles raises an FPE under -ftrapping-math at 4 MPI ranks
(clang+CUDA+FPE config) — reproduced identically at the pre-migration base
commit; upstream PR candidate alongside the particle_node prefetch predicate.

**Layout re-evaluation (2026-07-09, user-requested): COMPONENT-MAJOR becomes
the default — the phase-3.5 decision is REVERSED.** Method: a compile-time
toggle builds both layouts from identical sources (all column access was
already layout-agnostic except two hand-hoisted kernel pointer sites, made
stride-aware); trajectories are bitwise-identical under either layout, so the
comparison is purely speed. Order-balanced same-commit A/B (6 reps/config):
component-major is faster on p3m 4-rank (−1.7%), p3m 1-rank (−1.3%), lj
1-rank and 4-rank (−0.5% each); lj-omp was +1.0% and was recovered to +0.5%
(noise) by perf-guided optimization — an interleaved position scratch filled
once per Verlet build (the build's per-candidate reads touched three
far-apart component streams, ~2× its cost at 4000 ppc); p3m-omp is
inconclusive (readings 1.006 and 1.033 across rounds — P3M auto-tuning
variance dominates that configuration). Ship-gate adjudication: the
sequential gate read lj-1rank 1.055 vs phase-0 on decision day, but a
same-session order-balanced anchor measured the layout delta on lj-1rank as
exactly 1.000 (identical min-of-means) and BOTH layouts at 1.055 — i.e.
~0.6pp of cross-day machine drift affecting both binaries equally
(particle-major itself gated at 1.049 the previous day); per the established
A/B-arbitrates precedent the 5% budget verdict stands. Why the reversal: the
phase-3.5 measurement predated the struct death, per-field ghost exchange and
column kernels — the gather-dominated paths that favoured particle-major are
gone, and per-component ghost legs now benefit from contiguity.
Component-major is also the coalescing-friendly layout for future
device-resident GPU work. Particle-major remains selectable via
-DESPRESSO_PARTICLE_STORE_PARTICLE_MAJOR.

## Non-goals

- No change to the Python user interface or its semantics.
- No algorithmic changes to integrators, thermostats, or solvers — this is a
  storage migration; any numerical difference beyond reduction-order effects
  (phases 7–8) is a bug.
- Double precision remains canonical for all stored quantities; solvers that
  internally use single precision keep doing their own casting.
- MPI domain decomposition and the ghost-exchange communication topology are
  unchanged — only the pack/unpack representation changes.

## Section 1: Architecture & data model

**ParticleStore.** A new class (new directory `src/core/particle_store/`) owns
every per-particle quantity in a single index space: local particles first,
contiguously sorted by cell, ghost particles appended after `n_local`.
`Cell` degrades to an `(offset, count)` range into the store;
`Utils::Bag<Particle>` ownership disappears. A consequence of the layout:
`is_ghost(i)` becomes `i >= n_local` — the per-particle ghost flag is deleted.

```
index:   0 ................ n_local-1 | n_local .... n_total-1
order:   [cell 0][cell 1][cell 2]...  | [ghost cell 0][...]

Cell #k  = { offset, count }   // range into store
resort   = stable permutation of all columns
id2index = Kokkos::View<int*> (rebuilt on resort)
```

**Three field categories** (state/observables separation, refined one step):

1. **Parameters** — user-set, constant during integration: id, mol_id, type,
   propagation, rotation/fix bitfields, mass, charge, dipole moment,
   rotational inertia, gamma, gamma_rot, external force/torque, mu_E.
   Cold, feature-gated ones (swimming, Stoner-Wohlfarth, vs_relative) become
   host-only sidecar columns.
2. **State** — evolves during integration: position, image box, quaternion,
   velocity, omega, Lees-Edwards offset/flag, position at last Verlet update,
   RATTLE previous position.
3. **Observables** — recomputed or accumulated every step, never part of
   state proper: force, torque, RATTLE correction, dipole field. Written
   through `Kokkos::ScatterView` in parallel kernels, then contributed into
   the column.

Each hot field is a `FieldColumn<T, N>`: a DualView-style host/device pair of
`Kokkos::View<T*[N], LayoutLeft>` (component-major) plus per-field
modified/sync flags. In CPU-only builds host and device alias the same
memory — zero overhead. Scalars are `View<T*>`. Optional fields stay behind
the existing `ESPRESSO_*` ifdefs. Ragged data (bonds, exclusions) stays in
mutation-friendly host sidecars, flattened to Kokkos views on rebuild.

**One honest API cost:** with component-major layout there is no contiguous
`Vector3d` in memory, so today's reference-returning accessors (`p.pos()`)
cannot survive as-is. The proxy returns a small `VectorRef` supporting read,
`=`, `+=`, and conversion to `Utils::Vector3d` — most call sites compile
unchanged; the rest are mechanical fixes.

**id→index map:** `Kokkos::View<int*>` indexed by particle id (replacing
today's `std::vector<Particle*>` `m_particle_index`), rebuilt on resort;
serves Python access, bond partner resolution, and ghost packing.

## Section 2: Data flow

**Timestep loop** (all kernels are `Kokkos::parallel_for` over index ranges):

1. **Integrate** — velocity-Verlet kernel over `[0, n_local)`: reads
   parameters + observables (force), writes state (pos, v). Propagation-mode
   filtering via the propagation column, as today.
2. **Resort** (when triggered) — the decomposition computes each particle's
   target cell; a stable permutation is built and applied to every column;
   cell `(offset, count)` ranges, the id→index map, ghost rows, and flattened
   bond views are rebuilt. Resort already implies Verlet-rebuild today, so no
   new invalidation cadence is introduced.
3. **Ghost update** — per-field: gather send-rows into `CommBuf` by index,
   MPI exchange, scatter into ghost rows (indices ≥ `n_local`). The
   `GHOSTTRANS_*` flags map onto field groups: POSITION → pos+image+quat,
   MOMENTUM → v+omega, PROPRTS → parameter columns, FORCE → reverse path with
   reduction into local rows. Boost-serialization of `Particle` sub-structs is
   replaced by per-column memcpy-style gathers — bitwise, no archive overhead.
   Bonds-in-ghosts uses the sidecar as today.
4. **Neighbor lists** — Cabana `VerletList` builds directly from the position
   column. The entire `commit_particle` copy-in layer and the `AoSoA_pack`
   mirror are deleted — this is the payoff that funds the ≤5% regression
   budget.
5. **Forces** — zero/init observables columns (external forces, thermostat),
   short-range pair kernel via `Cabana::neighbor_parallel_for`, bonded kernels
   via the flattened bond views, long-range solvers (P3M) read charge/pos
   columns. All contributions go through ScatterView, then one contribute +
   ghost-force reduction.

**Python path** — unchanged UI. `ParticleHandle` resolves id→index on the
owning rank and reads/writes columns through the same MPI machinery as today.
Writes mark the affected field host-modified; the next device-side kernel
syncs only stale fields.

**GPU builds** — kernels launch in the default execution space; per-field
sync tracking means host-only subsystems (Python, ghosts before GPU-aware
MPI) pull only the fields they touch.

## Section 3: Migration phasing

The struct is not replaced in one step; it is **emptied field-by-field**.
Each phase moves one field group's source of truth from `Particle` into
`ParticleStore` columns; the `Particle` accessors for that group redirect
into the columns via a transitional `store_index` member, so all call sites
keep compiling and behavior is identical. When the struct is empty, it is
replaced by the view class. Every phase ends at a mergeable checkpoint:
full test suite green, LJ/P3M benchmarks within the 5% budget.

**Field-eviction order** (simplest ghost-communication direction first):

| Phase | Content | Checkpoint risk retired |
|---|---|---|
| 0 | Record LJ/P3M benchmark baselines for the agreed configs; add a comparison script; C++ unit-test scaffolding for `src/core/particle_store/`. | Perf gate exists before anything moves. |
| 1 | Preparatory refactor, no storage change: centralize every particle insert/remove/move behind a small set of `CellStructure` primitives (today `RegularDecomposition::resort` manipulates `ParticleList`s directly). | Single hook points for mirroring rows; the riskiest transitional machinery gets one home. |
| 2 | `ParticleStore` skeleton: row bookkeeping mirrors the phase-1 primitives; `store_index` member added to `Particle`; store order rebuilt (cell-sorted) at resort; topology changes mark the store dirty for rebuild. Evict **force + torque** (observables): columns become source of truth, `ScatterView` contribution, per-field ghost-force reduction. | Store/row-index consistency machinery proven on the field with the weakest persistence requirement (forces are recomputed after every resort, so resort may simply reallocate the column). |
| 3 | Evict **position, image box, quaternion** (+ `p_old`, `p_last_timestep`, Lees-Edwards offset/flag). Verlet/link-cell and short-range kernels read the position column directly; `commit_particle` position copies deleted; per-field POSITION ghost exchange. | The hot read path and the biggest copy-in cost. First phase where resort must permute state columns. |
| 4 | Evict **velocity + omega**. Thermostats and integrators via proxy accessors on columns; per-field MOMENTUM ghost exchange; remaining `commit_particle`/`AoSoA_pack` mirror deleted. | Dual-storage copy layer fully gone. |
| 5 | Evict **parameters** (id, mol_id, type, propagation, rotation/fix bitfields, mass, charge, dipm, rinertia, gamma, gamma_rot, ext_force/torque, mu_E) and move cold groups (swim, magnetodynamics, vs_relative) to host sidecar columns. Per-field PROPRTS ghost exchange; Python setters write columns. | id→index map replaces `m_particle_index`. |
| 6 | Evict **ragged data**: `BondList` and exclusions become sidecar columns indexed by store row; flattened kernel views unchanged; ghost-bond transfer from sidecar. | `Particle` is now empty except `store_index`. |
| 7 | **Kill the struct**: `Particle` becomes the view class (store pointer + index, same accessor API); `Utils::Bag<Particle>` in cells replaced by (offset, count) ranges; resort rewritten as pure column permutation; inter-rank particle exchange (global resort) switches from boost serialization to per-field pack; `GpuParticleData` replaced by device views of columns. | AoS storage no longer exists anywhere in the core. |
| 8 | **De-proxy hot paths & GPU enablement**: rewrite integrators, thermostats, and remaining per-particle loops as direct column kernels; validate CUDA build; enable device-resident execution with per-field sync. | End-state performance goals. |

**Transitional invariants** (checked in debug builds from phase 2 on):
`store.id(p.store_index()) == p.id()` for every particle after any resort or
topology change; column extents equal `n_local + n_ghost`; a field group is
owned by exactly one side (struct or column) at any commit — accessors of an
evicted field never touch struct memory.

Each phase is a separate sub-project with its own implementation plan
(phases 0 and 1 may share one); phase boundaries are the only merge points.

## Section 4: Error handling

**User-facing error behavior does not change.** Python-level exceptions
(accessing a nonexistent particle, invalid attribute values) and the runtime
error collection machinery (`errorhandling`) keep their exact semantics; bond
resolution/breakage errors already work from Kokkos kernels
(`bond_broken_error`) and are unaffected.

**Internal invariants become debug-build checks** rather than silent
assumptions:

- **Row consistency:** `store.id(p.store_index()) == p.id()` after every
  resort and topology change; column extents equal `n_local + n_ghost`.
  Adjudicated deviation (phase 2): the id-based cross-check
  (`store.id(p.store_index()) == p.id()`) is deferred until the id column
  exists (migration phase 5); phase 2 asserts only the row-index bounds
  (`0 <= row < n_local + n_ghost`) in `ParticleStore::assign_row`.
- **Single ownership:** an evicted field's struct accessor must never touch
  struct memory — enforced by removing the struct member in the same commit
  that evicts the field (compile-time, not runtime).
  Adjudicated deviation (phase 2, observables): migration-carrier members
  ferry force/torque values through whole-particle serialization (inter-rank
  migration, fetch cache); accessors never read them; `assign_row` seeds a
  detached arrival's row from its carriers.
  Adjudicated deviation (phase 3, state fields): for DETACHED particles
  (freshly constructed, deserialized, not yet assigned a row), state-field
  accessors read/write the migration carriers — detached state lives in
  carriers, attach seeds columns from carriers. Single ownership holds
  unconditionally for ATTACHED particles: their accessors touch only
  columns. New invariant this creates: never `assign_row` a particle into a
  foreign store without refreshing its carriers (serialize) or re-writing
  its state through accessors afterwards — stale carriers would seed the row.
  Adjudicated deviation (phase 3, kernels): `commit_particle`'s per-step
  position/image/director copies into the kernel pack were retained at the
  phase-3 flip (kernels still read the pack, not the store columns; the
  pack-index→store-row translation view was not built). Re-adjudication with
  benchmark data is owed at the phase-3 checkpoint; the copy layer is
  scheduled to die with the velocity eviction (phase 4) or a dedicated
  kernel-rewiring step.
- **Sync-state misuse:** in debug builds, reading a field in a memory space
  where it is stale aborts with the field name and the modifying site.
  In release builds sync is automatic (per-field, copy only when stale).
- **id→index lookups:** a missing id yields the sentinel `-1`, mapped to the
  same "particle not here" behavior as today's `nullptr` from
  `get_local_particle()`.
- **Ghost pack/unpack:** per-field byte counts are computed from the
  `GHOSTTRANS_*` flags on both sides; debug builds assert the received
  buffer size matches the expectation (replacing the implicit trust the
  boost-archive layer provided).
- **Allocation:** Kokkos host-side allocation failures propagate as
  exceptions before any column is partially resized (resize order:
  allocate-new, copy, swap).

## Section 5: Testing

- **Primary gate:** the full existing Python test suite (ctest) at every
  phase checkpoint, on 1 and 4 MPI ranks. No test may be weakened or skipped
  to make a checkpoint.
- **New C++ unit tests** for `src/core/particle_store/`: column permutation
  under resort, id→index rebuild, per-field DualView sync semantics,
  `VectorRef` proxy semantics (read/`=`/`+=`/conversion), per-field-group
  ghost pack/unpack roundtrips, and randomized add/remove/resort sequences
  validating the Section 3 invariants.
- **Sanitizers:** ASan/UBSan builds of the unit tests plus a representative
  subset of the Python suite once per phase (row-index bugs manifest as
  out-of-bounds column access — exactly what ASan catches).
- **Benchmarks:** LJ and P3M at `--particles_per_core 1000` (1 and 4 MPI
  ranks) and `--particles_per_core 4000` (4 OMP threads), compared against
  the phase-0 recorded baseline; budget ≤5% cumulative regression.
- **Shared-machine protocol:** the machine is shared. The benchmark
  comparison script checks load (`uptime`, foreign heavy processes) before
  measuring and refuses to record results on a loaded machine; runs are
  repeated several times and the minimum is used. No baseline or gate
  decision is taken from a busy machine.
- **Numerical identity check:** for phases 2–6 (pure storage moves, no
  algorithm change), a short LJ+P3M trajectory must be bitwise-identical
  (positions/forces) before and after the phase on a fixed seed, serial run.
  Phases 7–8 may relax this to statistical agreement where reduction order
  legitimately changes.
