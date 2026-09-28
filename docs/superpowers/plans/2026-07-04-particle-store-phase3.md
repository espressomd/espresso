# ParticleStore Migration — Phase 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evict position, image box, quaternion, `pos_at_last_verlet_update`, `pos_last_time_step`, and the Lees-Edwards offset/flag from the `Particle` struct into ParticleStore columns; short-range/bond kernels read the position column directly and the AoSoA position mirror dies.

**Architecture:** Proxies generalize (`BasicVectorReference<T>` for double/int 3-vectors, `QuaternionReference` for the 4-wide quaternion column; scalar columns return real element references). The kernel/layout risk is isolated from the ownership risk: Task 7 first swaps the kernels' position/image source from the AoSoA pack to store columns while the struct stays authoritative (a benchmarkable mirror swap), then Task 8 flips ownership (accessors → proxies, struct members deleted, per-step position commits deleted). Detached fetch-cache particles get a head-node snapshot store seeded from migration carriers (the adjudicated `part_rep` pattern), so cluster analysis and all ParticleHandle getters keep working unchanged. Ghost POSITION becomes value-based direction-first serialization (the phase-2 FORCE pattern) with compositional transmit sizing.

**Tech Stack:** C++20, Kokkos, Boost.MPI/Test, phase-0 benchmark gate, phase-2 identity gate.

**Spec:** `docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md` (sections 3 phase 3, 4, 5). **Phase-2 precedents referenced throughout:** migration carriers (`Particle.hpp` `detached_force()`/`m_detached_force`), value-preserving `assign_row`, sync points, `ParticleStoreTestFixture`, private stores for standalone particles (`ShapeBasedConstraint::m_part_rep_store`).

## Global Constraints

- Worktree `/tikhome/weeber/es/.claude/worktrees/eliminate_Particle` (branch `worktree-eliminate_Particle`) only; commands from worktree root; never touch other worktrees.
- `make -j8` only; no cmake re-configure; pre-commit hook may reformat (re-stage, re-commit); commit messages end with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Tests green at every commit. Phase gate: full Python suite; benchmarks ≤5% CUMULATIVE vs `maintainer/benchmarks/baselines/phase0-baseline.csv` via `benchmark_gate.py` (foreign-CPU busy check; shared machine — never bypass).
- **Identity gate at every task that touches numerics paths:** both modes of `maintainer/benchmarks/trajectory_identity.py` vs `.superpowers/sdd/phase1-identity-reference.txt` (lj `bdd2022c…`, p3m `141cc4aa…`) must be IDENTICAL. A mismatch is BLOCKED, never rationalized.
- **Single ownership (spec section 4):** the flip commit (Task 8) deletes the struct members in the same commit that makes columns authoritative. Carriers are serialization envelopes only — accessors never read them.
- Test gate patterns MUST include `lees_edwards|hybrid_decomposition|virtual_sites` for any storage/ghost change (phase-2 lesson).
- After `Particle.hpp` ABI changes: if inexplicable segfaults/zero-physics appear, `make -C build clean && make -C build -j8` FIRST (NFS clock-skew stale objects; phase-2 lesson).
- Known ignorable test noise: `sample_*`, `tutorial_*` (separate fixture suites), `benchmark_lb_with_cpu`. Anything else is real.
- Full-word naming. GPL headers on new files.
- HDF5 side-build check at the checkpoint (phase 2 broke HDF5 linking once): configure a lean HDF5-on build only if `h5md_core.cpp`-adjacent files were touched — Task 9 decides.

---

### Task 1: Generalize proxies — `BasicVectorReference<T>` and `QuaternionReference`

**Files:**
- Modify: `src/core/particle_store/VectorReference.hpp`
- Test: `src/core/unit_tests/VectorReference_test.cpp` (extend)

**Interfaces:**
- Produces (binding for Tasks 2, 7, 8):
  - `template <class T> class BasicVectorReference` — same API as today's `VectorReference` but element type `T`: ctor `(T *base, std::size_t stride)`, `operator Utils::Vector<T,3>() const`, `operator=(Utils::Vector<T,3> const&)`, `+=`, `-=`, `*=(T)`, `operator[](std::size_t)->T&` + const, `norm()`, `norm2()` (arithmetic members only when `T` is floating point is NOT required — keep them; they compile for int too), write-through copy assignment (phase-2 convention).
  - `using VectorReference = BasicVectorReference<double>;` and `using IntegerVectorReference = BasicVectorReference<int>;`
  - `class QuaternionReference` — over a 4-wide double column: ctor `(double *base, std::size_t stride)`, `operator Utils::Quaternion<double>() const`, `operator=(Utils::Quaternion<double> const&)`, `operator+=(Utils::Quaternion<double> const&)`, `operator/=(double)`, `operator[](std::size_t)->double&` + const (index 0..3), `norm()`, `norm2()`, write-through copy assignment.

- [ ] **Step 1: Write failing tests** — extend `VectorReference_test.cpp` with: an `IntegerVectorReference` write/read/subscript case over a stride-4 int column fixture; a `QuaternionReference` case (stride-4 double column of 4 components: assignment from `Utils::Quaternion<double>`, `operator[]` component write, `+=`, `/=`, conversion round-trip, write-through copy assignment between two rows).
- [ ] **Step 2: Run to verify failure** (`make -C build -j8 VectorReference_test` → compile errors).
- [ ] **Step 3: Implement** — templatize the existing class body over `T` (mechanical: `double`→`T`, `Utils::Vector3d`→`Utils::Vector<T,3>`); keep the class doc (template-argument-deduction warning). Add `QuaternionReference` as a separate class in the same header (4 components; `Utils::Quaternion<double>` exposes `operator[]`/`.norm()` — check `src/utils/include/utils/quaternion.hpp` for the exact element access API and mirror it).
- [ ] **Step 4: Tests pass** (`ctest --test-dir build -R VectorReference_test`).
- [ ] **Step 5: Commit** (`core: generalize store proxies to BasicVectorReference<T> and QuaternionReference (phase 3)`).

---

### Task 2: ParticleStore state columns + generalized rebuild machinery

**Files:**
- Modify: `src/core/particle_store/ParticleStore.hpp`, `src/core/particle_store/ParticleStore.cpp`
- Modify: `src/core/Particle.hpp` (state migration carriers)
- Test: `src/core/unit_tests/ParticleStore_test.cpp` (extend)

**Interfaces:**
- Produces (binding for Tasks 7, 8):
  - New columns + accessors on `ParticleStore`:
    - `position_reference(int)->VectorReference` / `position_value(int)->Utils::Vector3d`; column type `Kokkos::DualView<double*[3], Kokkos::LayoutLeft>` named `particle_store::position`
    - `image_box_reference(int)->IntegerVectorReference` / `image_box_value(int)->Utils::Vector3i`; `DualView<int*[3], LayoutLeft>`
    - `quaternion_reference(int)->QuaternionReference` / `quaternion_value(int)->Utils::Quaternion<double>`; `DualView<double*[4], LayoutLeft>` (ESPRESSO_ROTATION)
    - `position_at_last_verlet_update_reference(int)->VectorReference` / `..._value(int)`
    - `position_last_time_step_reference(int)` / `..._value(int)` (ESPRESSO_BOND_CONSTRAINT)
    - `lees_edwards_offset(int)->double&` and `lees_edwards_flag(int)->short&` (scalar columns `DualView<double*>` / `DualView<short*>` — host element reference directly, no proxy)
    - `position_view()` etc. host-view getters for kernel wiring (Task 7).
  - Refactor `begin_rebuild`/`assign_row`/`finish_rebuild` to iterate a per-column generalized copy (avoid 9× copy-paste): a private helper template applied to each column pair, e.g. `preserve_or_seed(new_column, old_column, row, old_row, carrier_value)`. Signatures of the three public functions unchanged.
  - Carriers on `Particle` following the exact phase-2 pattern (`detached_force`): `detached_position()`, `detached_image_box()`, `detached_quaternion()`, `detached_position_at_last_verlet_update()`, `detached_position_last_time_step()`, `detached_lees_edwards_offset()`, `detached_lees_edwards_flag()` + `m_migration_*` members, serialized in `Particle::serialize` next to the force carriers; `assign_row` seeds all state columns from carriers for detached particles. NOTE: carriers are IN ADDITION to the still-present struct members until Task 8 — during Tasks 2-7 the carrier getters read the struct fields (attached or not), and `Particle::serialize` keeps serializing the ORIGINAL sub-structs; the carrier members/serialization entries are added but DORMANT until Task 8 rewires them (add them now so the serialization format changes only once, in Task 8: therefore do NOT add carrier fields to `serialize` in this task — only the members and getters; Task 8 adds the `ar &` lines when the sub-structs leave).
- [ ] **Step 1: Failing tests** — extend `ParticleStore_test.cpp`: rebuild preserves position/image/quaternion/le-offset by old row (set via references, shuffle rows as in the phase-2 test, assert values followed the particles); detached particle with carrier members set gets its state columns seeded on `assign_row` (construct a `Particle`, set its migration carrier members via a test hook or friend — check how the phase-2 test `rebuild_seeds_migrated_particle_force_from_carrier` sets carriers and mirror it); scalar column references write through.
- [ ] **Step 2: verify failure.**
- [ ] **Step 3: Implement** per the Interfaces block. Bounds assert in `assign_row` already exists (phase-2 fix) — keep.
- [ ] **Step 4: All unit tests green** (`ctest --test-dir build -L unit_test`).
- [ ] **Step 5: Commit** (`core: ParticleStore state columns and generalized rebuild (phase 3)`).

---

### Task 3: Call-site preparation audit A — positions, images, fold sites, ranges

**Files:** (grep-driven; known criticals listed)
- Modify: `src/core/cell_system/RegularDecomposition.cpp` (`fold_and_reset`), `src/core/cell_system/AtomDecomposition.cpp:110-112`, `src/core/lees_edwards/lees_edwards.hpp:69`, `src/core/virtual_sites/relative.cpp:147`
- Modify: `src/core/ParticlePropertyIterator.hpp:54` (+ its consumer `src/core/electrostatics/elc.cpp:1216-1258`)
- Modify: every site the greps below surface.

**Interfaces:** zero behavior change; today's reference-returning accessors keep working. Produces a codebase where `pos()`, `image_box()`, `pos_at_last_verlet_update()`, `pos_last_time_step()`, `lees_edwards_*()` appear only in proxy-compatible patterns.

- [ ] **Step 1: fold_position sites** — all 4 known callers (plus any more from `grep -rn "fold_position(p\.\|fold_position(part\." src/core`):
```cpp
    Utils::Vector3d position = p.pos();
    Utils::Vector3i image_box = p.image_box();
    box_geo.fold_position(position, image_box);
    p.pos() = position;
    p.image_box() = image_box;
```
  (`fold_and_reset` additionally keeps `p.pos_at_last_verlet_update() = p.pos();` — reads the fresh value, still fine.)
- [ ] **Step 2: pos_range rework** — `ParticlePropertyIterator.hpp`: change `return_pos` lambda to return `Utils::Vector3d` BY VALUE for now (`[](Particle &p) -> Utils::Vector3d { return p.pos(); }`); check elc.cpp consumers (`charge_assign`, `modify_p3m_sums`, and the p3m `charge_assign` signature —
  `grep -n "charge_assign\|prepare_sc_cache" src/core/electrostatics/elc.cpp src/core/electrostatics/p3m*.hpp`) accept values/const& (they read positions; confirm no consumer WRITES through pos_range — if one does, report it in the task report as a deviation and fix write sites individually). NOTE: post-flip (Task 8) this lambda becomes `-> Utils::Vector3d { return p.pos(); }` reading the proxy — the value-return form chosen here survives the flip unchanged.
- [ ] **Step 3: greps for remaining hostile patterns** (fix every hit with explicit `Utils::Vector3d(...)` copies / explicit-type locals; the phase-2 Task-4 report shows the pattern catalog):
```bash
grep -rn "auto &[a-z_]* = [a-zA-Z0-9_]*\.pos()\|auto &[a-z_]* = [a-zA-Z0-9_]*\.image_box()" src/core src/script_interface
grep -rn "\.pos() [-+*/]" src/core src/script_interface | grep -v "+=\|-=\|\*=\|/="
grep -rn "auto [a-z_]* = [a-zA-Z0-9_]*\.pos()\b" src/core src/script_interface
grep -rn "\.pos_at_last_verlet_update() [-+]" src/core | grep -v "="
```
  Known template-arithmetic sites needing explicit conversion: `CellStructure.cpp:664` verlet criterion `(p.pos() - p.pos_at_last_verlet_update()).norm2()`; `lb_tracers.cpp:93` same shape; statistics/analysis subtraction chains. Also `System.cpp:439` component access.
- [ ] **Step 4: Gates** — build; unit tests; `ctest -R "elc|coulomb|lees_edwards|cellsystem|^particle$|analysis|cluster"`; identity BOTH modes IDENTICAL.
- [ ] **Step 5: Commit** (`core: prepare position/image call sites for proxy accessors (phase 3)`).

---

### Task 4: Call-site preparation audit B — quaternions, director, dipoles

**Files:** `src/core/rotation.cpp`, `src/core/rotation.hpp`, `src/core/virtual_sites/relative.cpp`, `src/script_interface/particle_data/ParticleHandle.cpp`, plus grep hits.

**Interfaces:** zero behavior change; `quat()` sites reduced to: `= q`, `+= q` (explicit-converted where template deduction fails), `[i]` component access, passed as `Quaternion const&` (conversion), `.norm()/.norm2()`.

- [ ] **Step 1: greps**
```bash
grep -rn "auto &[a-z_]* = [a-zA-Z0-9_]*\.quat()" src/core src/script_interface
grep -rn "\.quat() [-+*/]" src/core src/script_interface | grep -v "+=\|/="
grep -rn "\.quat()\." src/core src/script_interface
```
- [ ] **Step 2: Known sites** — `rotation.cpp:152` `p.quat() += time_step * (Qd + ...) - lambda * p.quat()`: RHS uses `p.quat()` in template arithmetic → `auto const quaternion_old = Utils::Quaternion<double>(p.quat());` local, compute RHS from `quaternion_old`, then `p.quat() += rhs;` — verify the existing math reads the PRE-update quat consistently (it does — single statement). `rotation.cpp:159` `p.quat() /= scale` — proxy has `/=`, fine, leave. Component reads `p.quat()[0..3]` in `define_Qdd` (rotation.cpp:65-110) take `Particle const&`? Check — if const context, post-flip const accessor returns a VALUE `Quaternion` — component reads fine; if it's called per-component many times in a hot loop, hoist one `auto const quaternion = ...` local (do it now for clarity). `relative.cpp:153` `p.quat() = p_ref.quat() * p.vs_relative().quat` — RHS product of proxy: `Utils::Quaternion<double>(p_ref.quat()) * ...` explicit.
- [ ] **Step 3: `convert_quaternion_to_director(p.quat())` sites** (commit_particle, calc_director) — takes `Quaternion const&` → conversion handles; verify signature.
- [ ] **Step 4: Gates** — build; unit tests; `ctest -R "rotation|dipol|magneto|virtual_sites|engine"`; identity IDENTICAL.
- [ ] **Step 5: Commit** (`core: prepare quaternion call sites for proxy accessors (phase 3)`).

---

### Task 5: Ghost POSITION value-based serialization

**Files:** `src/core/ghosts.cpp`

**Interfaces:** wire format unchanged (pos 3d + image 3i + quat 4d [ROTATION] + p_last_time_step 3d [BOND_CONSTRAINT], in that order). Produces a POSITION branch that never binds particle struct memory to archives; compositional transmit size for POSITION (mask-out like FORCE).

- [ ] **Step 1: Rework `serialize_and_reduce` POSITION branch** direction-first (phase-2 FORCE pattern, but position has NO reduction policy — MOVE semantics both policies):
```cpp
  if (data_parts & GHOSTTRANS_POSITION) {
    if (direction == SerializationDirection::LOAD) {
      Utils::Vector3d position;
      Utils::Vector3i image_box;
      ar & position;
      ar & image_box;
      p.pos() = position;
      p.image_box() = image_box;
    } else if (ghost_shift != nullptr) {
      /* ghost shift + fold on SAVE (existing semantics, already value-based) */
      Utils::Vector3d position = p.pos() + *ghost_shift;
      Utils::Vector3i image_box = p.image_box();
      box_geo.fold_position(position, image_box);
      ar & position;
      ar & image_box;
    } else {
      Utils::Vector3d position = p.pos();
      Utils::Vector3i image_box = p.image_box();
      ar & position;
      ar & image_box;
    }
#ifdef ESPRESSO_ROTATION
    /* same direction-first value pattern for p.quat() (4 doubles) */
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
    /* same for p.pos_last_time_step() */
#endif
  }
```
  READ THE CURRENT BRANCH FIRST (ghosts.cpp:172-190) and preserve every semantic including the shift+fold path exactly. Lees-Edwards: verify by grep whether the POSITION branch or the shift computation touches `lees_edwards_offset/flag` (`grep -n "lees_edwards" src/core/ghosts.cpp src/core/cell_system/RegularDecomposition.cpp`) — if the LE offset/flag are ghost-communicated anywhere, handle them the same value-based way and note it in the report.
- [ ] **Step 2: `calc_transmit_size`** — add POSITION to the mask-out scheme: compute `position_size` as a compile-time expression (3*sizeof(double) + 3*sizeof(int) + [4*sizeof(double)] + [3*sizeof(double)]) — CAREFUL: verify against the CURRENT sizer output first (`MemcpyOArchive` may pad `Vector3i` — write a quick throwaway check or read `Utils::MemcpyOArchive` to confirm packing is tight); mask POSITION out of `data_parts` before the `Particle{}` sizer, add `position_size` to the return. Keep FORCE handling as-is.
- [ ] **Step 3: Gates** — build; unit tests; `ctest -R "particle|cellsystem|coulomb|lees_edwards|hybrid|virtual_sites|npt"`; identity IDENTICAL (ghost positions feed forces directly — any wire slip shows).
- [ ] **Step 4: Commit** (`core: ghost position exchange via explicit field values (phase 3)`).

---

### Task 6: Fetch-cache snapshot store + sync choke point

**Files:** `src/core/particle_node.cpp` (+ `.hpp` if a helper needs declaring)

**Interfaces:**
- Produces: every `Particle` placed into `particle_fetch_cache` is attached to a head-node-local `ParticleStore` ("snapshot store") seeded from its migration carriers, so ALL accessor reads on cached particles keep working after the flip (cluster analysis, ParticleHandle getters). Pattern precedent: `ShapeBasedConstraint::m_part_rep_store` + `assign_row` carrier seeding.
- Mechanism:
  - File-scope `static ParticleStore fetch_cache_store;` + row counter next to `particle_fetch_cache`.
  - A helper `attach_cached_particle(Particle &p)` that: if the store is full or unbuilt, grows it (`begin_rebuild(capacity, 0)` / `finish_rebuild()` with doubling capacity — values need not survive growth EXCEPT rows already attached must be preserved; simplest correct scheme: on grow, rebuild with the doubled capacity and re-`assign_row` every currently cached particle in cache order — enumerate via a side vector of cached pids/particle pointers maintained alongside; keep it simple and correct, this is a cold path); then `assign_row(p, next_row++)` (seeds from carriers).
  - Call it on BOTH cache-fill paths: `get_particle_data`'s `particle_fetch_cache.put(...)` (the returned pointer's particle) and the prefetch bulk path (`fetch_particles`' puts, ~line 279).
  - `invalidate_fetch_cache()` additionally resets the snapshot store and row counter (cache and store lifecycles stay in lock-step).
  - Sync choke point: first statement of `mpi_send_particle_data_local` AND of `get_particle_data`'s local-return path: `get_cell_structure().ensure_particle_store_synchronized();` (owner-side serialize reads accessors for carrier fill → owner store must be clean; the local-return path hands out a live particle whose reads need valid rows).
- Pre-flip (now), attachment is harmless (accessors still read structs) — this task lands the machinery and its tests while everything is still observable.
- [ ] **Step 1: Implement** as above.
- [ ] **Step 2: Gates** — build; `ctest -R "cluster|analysis|^particle$|particle_slice"` (cluster analysis exercises cached reads across ranks — run the 4-rank variants); unit tests; identity IDENTICAL.
- [ ] **Step 3: Commit** (`core: fetch-cache snapshot store and owner-side sync choke point (phase 3)`).

---

### Task 7: Kernel mirror swap — kernels read store position/image columns

**Files:** `src/core/aosoa_pack.hpp`, `src/core/short_range_cabana.hpp`, `src/core/forces_cabana.hpp`, `src/core/energy_cabana.hpp`, `src/core/pressure_cabana.hpp`, `src/core/bond_forces_kokkos.hpp`, `src/core/bond_energy_kokkos.hpp`, `src/core/bond_pressure_kokkos.hpp`, `src/core/cell_system/CellStructure.hpp` (kernel-facing store view plumbing), `src/core/system/System.cpp` (`rebuild_aosoa` callers if signatures shift)

**Interfaces:**
- The struct stays authoritative (flip is Task 8); `commit_particle` now writes position/image into the STORE columns (via `store.position_reference(row) = p.pos()` or, better, direct host-view writes indexed by the same aosoa index — READ CAREFULLY: the aosoa index (`m_unique_particles` order) is NOT the store row! `commit_particle` receives the pack index; the store columns are indexed by store row. **Resolution: the kernels must index positions by the SAME index as forces (`m_local_force` uses the aosoa/unique-particles index). Therefore this task gives the PACK new position/image views that alias… NO — decision: keep it simple and correct:** the pack keeps its own position/image Kokkos Views BUT their layout changes to `LayoutLeft` and their element type/shape matches the future store columns exactly; `commit_particle` keeps writing them per-step from `p.pos()`. Kernels are updated where component-order/layout assumptions leak (audit `aosoa.position(i, j)` sites — Kokkos View indexing is layout-agnostic at the source level, so most sites compile unchanged; `get_span_at` (contiguous span over one particle's row) is LayoutRight-dependent — find its users (`grep -n get_span_at src/core -r`) and convert them to `get_vector_at`).
- So Task 7 = LayoutLeft conversion of the pack's position/image/director/velocity views + `get_span_at` elimination + benchmark spot-check. The actual "kernels read store columns" unification happens in Task 8 when the store columns become the live source and `commit_particle`'s position/image/director copies are DELETED with kernels pointed at store views by aosoa-index→store-row translation... **STOP — this needs the index-mapping resolution:**
- **Index mapping resolution (binding for Task 8):** kernels index particles by pack index `i` (from `m_unique_particles`); the store indexes by row. After the flip, kernels must read positions via `store_row_of_pack_index` — an int view built during `set_index_map` (`store_row_of[i] = m_unique_particles[i]->store_row()`), OR the verlet list build switches to store rows entirely. **Decision for this plan: build the translation view** (`Kokkos::View<int*> m_pack_index_to_store_row`, filled in `set_index_map` after `ensure_particle_store_synchronized()` — one extra indirection per kernel particle access for positions, benchmark-gated). If the benchmark gate fails on this, escalate per BLOCKED protocol (the alternative — reordering the store to unique_particles order — is a design change needing review).
- [ ] **Step 1: LayoutLeft conversion of pack views** (`aosoa_pack.hpp`: `PositionViewType` etc. → `Kokkos::LayoutLeft`), delete `get_span_at` after converting its users to `get_vector_at`.
- [ ] **Step 2: Build + full kernel test battery** — `ctest -R "coulomb|dipol|lj|gay|bond|angle|dihedral|interactions|nonbonded|pressure|energy|virtual_sites|lees_edwards|hybrid"`; unit tests; identity IDENTICAL.
- [ ] **Step 3: Benchmark spot-check (quiet machine only):** `benchmark_gate.py check-load` then `run --repetitions 2 --output /tmp/phase3-mirror.csv` and `compare` vs phase-0 baseline. Layout change is THE perf risk of this phase — if >5%, STOP and report BLOCKED with the table (do not proceed to the flip on a losing layout).
- [ ] **Step 4: Commit** (`core: component-major layout for kernel particle views (phase 3)`).

---

### Task 8: THE FLIP — position group moves to ParticleStore columns

**Files:** `src/core/Particle.hpp` (delete `ParticlePosition r` members pos/image/quat/p_last_timestep, `ParticleLocal` p_old + lees_edwards fields; accessors → proxies; serialize rewired to carriers), `src/core/cell_system/CellStructure.{hpp,cpp}` (`set_index_map` builds `m_pack_index_to_store_row`; kernel plumbing hands store position/image views + translation view to kernels), `src/core/short_range_cabana.hpp` (delete position/image commits from `commit_particle`; director commit now reads quat proxy), all kernel headers (position reads via translation view), `src/core/ghosts.cpp` (nothing — Task 5 prepared it), unit tests (fixture adaptations), plus the compile sweep.

**Interfaces:**
- `Particle::pos()` non-const → `VectorReference`; const → `Utils::Vector3d` value. `image_box()` → `IntegerVectorReference` / `Utils::Vector3i`. `quat()` → `QuaternionReference` / value. `pos_at_last_verlet_update()`, `pos_last_time_step()` → `VectorReference`/value. `lees_edwards_offset()` → `double&` (store scalar column element), const → `double`. `lees_edwards_flag()` → `short&`/`short`. All with `assert(m_particle_store != nullptr)`.
- `Particle::serialize`: `ar & r;` and the `ParticleLocal` position fields leave; carriers added (`ar & m_migration_position; ...` — SAVE fills carriers from `detached_*()` getters first, exactly like the force carriers at Particle.hpp:725).
- `calc_director()` / `calc_dip()` keep working (read quat by value).
- Ghost-cell particles, `part_rep`, fetch-cache particles: already covered (stores exist).
- Sync points: audit — position reads happen in resort itself (`particle_to_cell(p)` reads pos DURING resort while rows are valid — OK: dirty is marked after; verify no position access between extract-buffer-insert and the sync (buffered particles are detached only across MPI — local moves keep rows). New sync sites needed: `System::on_observable_calc` if not already present (check), `lees_edwards` protocol setter paths, `analysis` entries reading positions on LIVE particles outside PidObservable (grep callers of `local_particles()` in src/core/analysis/): add `ensure_particle_store_synchronized()` where missing.
- Unit-test fixtures: extend `ParticleStoreTestFixture` usage to every test constructing Particles that touch pos/quat/image (compile+runtime sweep will enumerate; known: rotation_test, Verlet_list_test, lees_edwards_test, link_cell_test, ParticleIterator/Bag tests if they set positions, BoxGeometry tests using Particle? — sweep decides).
- [ ] **Step 1: Flip Particle.hpp** (members out, proxies in, serialize rewired with carriers).
- [ ] **Step 2: Compile sweep** — iterate `make -C build -j8 2>&1 | grep error | head -40`, fixing per Task-3/4 rules; wire `m_pack_index_to_store_row` + kernel position/image view rewiring; DELETE position/image (+director-input change: director computed from `Utils::Quaternion<double>(p.quat())`) from `commit_particle`; delete the pack's position/image views.
- [ ] **Step 3: Runtime sweep** — unit tests (fixture adaptations); then the FULL pattern battery: `ctest -R "particle|cellsystem|integrator|langevin|brownian|npt|rotation|coulomb|dipol|lees_edwards|virtual_sites|collision|lb|cluster|analysis|hybrid|rattle|constraints|engine|dpd|observable|accumulator"`.
- [ ] **Step 4: Identity gate** — BOTH modes IDENTICAL (positions are hashed directly; this is the strongest signal of the phase).
- [ ] **Step 5: Commit** (`core: position, image box, and quaternion live in ParticleStore columns (phase 3 flip)`).

---

### Task 9: Phase-3 checkpoint

- [ ] **Step 1:** Guard script; unit tests; `hybrid_decomposition` repeat-until-fail:3.
- [ ] **Step 2:** `make -C build -j8 check_python_skip_long` then `check_python` — fully green (retry-once rule for statistical tests; deterministic failure = BLOCKED).
- [ ] **Step 3:** Identity both modes IDENTICAL.
- [ ] **Step 4:** Benchmark gate on final HEAD (quiet machine, 3 repetitions) vs phase-0 baseline — PASS ≤5% cumulative; report full table + call out the ratios for lj (position path is the hot path of this phase).
- [ ] **Step 5:** HDF5 side-build check: `h5md` files were not touched by this phase unless the sweep touched them — if `git diff bb9072bf..HEAD --stat | grep -E "h5md|io/writer"` shows hits, run the lean HDF5-on verification build from the phase-2 fix wave's recipe (see `.superpowers/sdd/phase2-final-fixes-report.md`).
- [ ] **Step 6:** Ledger append; nothing to commit.

---

## Plan Self-Review Notes

- **Spec coverage:** position/image/quat + p_old + p_last_timestep + LE fields evicted (Tasks 2, 8); kernels read the position column directly with the pack's position mirror deleted (Tasks 7-8 — via the pack-index translation view; the `commit_particle` position copies die in Task 8); per-field POSITION ghost exchange incl. shift + LE audit (Task 5); resort-must-preserve-state (value-preserving rebuild, Task 2; carriers for cross-rank, Tasks 2, 8); sync choke-point consolidation (Task 6 + Task 8 audit); identity + benchmark gates (global constraints, Task 9).
- **Known open risk, escalation defined:** kernel position reads gain one `pack-index → store-row` indirection (Task 7 Interfaces); benchmark-gated with BLOCKED escalation. The alternative (store ordered by unique_particles) is a design change → user input.
- **Type consistency:** `BasicVectorReference<T>`/`IntegerVectorReference`/`QuaternionReference` (T1) used in T2 accessors and T8 flip; `m_pack_index_to_store_row` named identically in T7/T8; carrier getters `detached_*` mirror phase-2 naming.
- Line numbers as of `8c852f88af`; code snippets authoritative over line numbers.
