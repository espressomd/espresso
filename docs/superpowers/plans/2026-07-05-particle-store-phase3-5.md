# ParticleStore Migration — Phase 3.5 Implementation Plan (kernel/store unification)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Kernels and the Verlet-list build read the ParticleStore position/image columns directly; `commit_particle`'s position/image/director copies are deleted; hot loops stop re-materializing proxies — reclaiming the phase-3 benchmark regression to pass the ≤5% gate so phase 3 can ship.

**Architecture:** Pack index and store row provably coincide for LOCAL particles (both assigned in cell-traversal order — `enumerate_local_particles` exclusive-scan vs `ensure_particle_store_synchronized`; assert it). Only the deduplicated ghost tail differs, so a `Kokkos::View<int*> m_pack_index_to_store_row` (identity for `i < n_local`, per-particle `store_row()` for the ghost tail) built in `set_index_map` translates kernel indices. The pack's position/image/director views die; kernels take the store's host views + the translation view. Director becomes a store-side derived view recomputed from the quaternion column in one bulk `parallel_for` per commit step (Gay-Berne/dipoles only). Integrator/thermostat hot loops hoist one proxy per field per particle instead of per component.

**Tech Stack:** C++20, Kokkos, phase-0 benchmark gate, phase-2 identity gate.

**User decision context (2026-07-04):** phase-3 checkpoint failed the benchmark gate (lj-4rank 1.078, p3m-1rank 1.059, p3m-4rank 1.077); user chose kernel unification over budget relaxation. Phase 3 is NOT pushed until this phase's gate passes.

## Global Constraints

- Same as phase 3 (worktree-only, `make -j8`, identity gate at every task — lj `bdd2022c…` / p3m `141cc4aa…` from `.superpowers/sdd/phase1-identity-reference.txt`, clean-rebuild-after-ABI-change rule, test noise allowlist, full-word naming, commit trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`).
- Benchmarks only on a quiet machine via `benchmark_gate.py`; ALL benchmark invocations run from the CONTROLLER session, never from subagents (background processes die with a subagent's turn — three incidents).
- Kernel index-space semantics must not change: ghost dedup (`m_unique_particles`), minimum-image convention, ScatterView force accumulation all stay as they are.

---

### Task 1: Kernels read store columns; commit_particle position copies die

**Files:**
- Modify: `src/core/cell_system/CellStructure.{hpp,cpp}` (`set_index_map` builds the translation view + local-prefix assertion; expose store position/image host views + translation view + director view to kernel plumbing)
- Modify: `src/core/aosoa_pack.hpp` (delete position/image/director views; keep velocity + scalar views)
- Modify: `src/core/short_range_cabana.hpp` (commit_particle loses position/image/director writes; verlet build already row-based from T8 — verify; wire new views into kernel data structs)
- Modify: `src/core/forces_cabana.hpp`, `src/core/energy_cabana.hpp`, `src/core/pressure_cabana.hpp`, `src/core/bond_forces_kokkos.hpp`, `src/core/bond_energy_kokkos.hpp`, `src/core/bond_pressure_kokkos.hpp` (~50 position reads + 1 image read + 4 director reads: `aosoa.position(i,j)` → `position(row(i),j)`; `get_vector_at(aosoa.position, i)` → `get_vector_at(position, row(i))` — keep `get_vector_at` as a free/static helper or member on a slim pack)
- Modify: `src/core/system/System.cpp` / ICC path if it touches the pack position views (grep `rebuild_aosoa` consumers)

**Interfaces:**
- Produces: `CellStructure::pack_index_to_store_row()` (`Kokkos::View<int const*>`), `store_position_view()`, `store_image_view()` (host views of the DualView columns), `director_view()` (derived, `View<double*[3], LayoutLeft>`, sized `n_total`, recomputed by `update_director_view()` — a `parallel_for` over rows converting the quaternion column; called where commit_particle's director write used to happen, guarded by the same GB/dipoles ifdefs).
- Invariant (assert in debug at build time of the translation view): `pack_index_to_store_row(i) == i` for all `i < count_local_particles()`.

- [ ] **Step 1:** Build translation view in `set_index_map` (after `ensure_particle_store_synchronized()` — add that call if absent): loop `m_unique_particles`, `view(i) = m_unique_particles[i]->store_row()`; debug-assert identity prefix.
- [ ] **Step 2:** Director derived view + `update_director_view()`; delete director from pack & commit.
- [ ] **Step 3:** Rewire kernel headers (mechanical; the translation lookup `row(i)` hoisted ONCE per particle index per kernel body — never per component).
- [ ] **Step 4:** Delete position/image from pack & commit_particle; fix `bond_forces_kokkos.hpp:239` image read via store image view + translation.
- [ ] **Step 5:** Clean core rebuild (ABI-adjacent header churn: `rm -rf build/src/core/CMakeFiles/espresso_core.dir`, full `make -C build -j8`, `make -C build -j8 unit_tests_executables`).
- [ ] **Step 6:** Gates: unit tests ALL; battery `ctest -R "particle|cell_system|coulomb|dipol|gay|nonbonded|interactions|bond|angle|dihedral|pressure|energy|virtual|lees_edwards|hybrid|collision|dpd|npt|rattle|lb" --output-on-failure`; identity BOTH modes IDENTICAL.
- [ ] **Step 7:** Commit `core: kernels read ParticleStore columns directly (phase 3.5)`.
- [ ] **Step 8:** CONTROLLER runs the quick benchmark (2 reps) and decides continuation.

### Task 2: Hot-loop proxy hoisting

**Files:** `src/core/integrators/velocity_verlet_inline.hpp`, `symplectic_euler_inline.hpp`, `brownian_inline.hpp`, `velocity_verlet_npt_*.cpp`, `steepest_descent.cpp`, `src/core/thermostats/langevin_inline.hpp`, `src/core/lb/particle_coupling.cpp`, `src/core/virtual_sites/{relative,lb_tracers}.cpp`, `src/core/rotation.cpp` (already partly hoisted), plus any `p.pos()[j]` / `p.v()[j]` / `p.force()[j]` per-component-in-loop patterns found by grep in hot inline headers.

**Interfaces:** none new — mechanical hoists: ONE proxy (or one value read + one write-back) per field per particle; component loops index the hoisted proxy/local. Semantics identical; identity gate must stay bitwise.

- [ ] **Step 1:** Grep-driven hoist pass over the listed files: pattern `for (j) { ... p.field()[j] ... }` → `auto field = p.field();` before the loop (proxy construction hoisted; `field[j]` hits memory directly). Where a field is only READ in arithmetic, hoist a `Utils::Vector3d const` value instead.
- [ ] **Step 2:** Gates: unit tests; `ctest -R "integrator|langevin|brownian|npt|steepest|lb|virtual|rotation"`; identity IDENTICAL.
- [ ] **Step 3:** Commit `core: hoist per-particle proxies in hot loops (phase 3.5)`.
- [ ] **Step 4:** CONTROLLER quick benchmark (2 reps): if ≤5% on all six, proceed to Task 3; if not, controller decides (profile deeper or escalate).

### Task 3: Phase-3(+3.5) checkpoint and ship

- [ ] **Step 1:** Guard; unit tests; hybrid repeat-until-fail:3.
- [ ] **Step 2:** Full `check_python_skip_long` + `check_python` — green.
- [ ] **Step 3:** Identity both modes IDENTICAL.
- [ ] **Step 4:** CONTROLLER: benchmark gate, 3 repetitions, vs phase-0 baseline — MUST PASS ≤5% on all six configs.
- [ ] **Step 5:** Update spec's phase-3 adjudication note (kernel rewiring now done; commit_particle keeps only velocity). Ledger. Push phase 3 + 3.5 to `me`.

## Self-Review Notes
- The ghost-tail translation gather is the only new per-access cost; locals (the bulk) hit the identity prefix. If Task-1 benchmarks show the gather hurting, the fallback (reorder m_unique_particles ghost tail to store-row order — dedup semantics preserved since dedup selects WHICH instances participate, not their order) is a contained follow-up.
- Type consistency: `pack_index_to_store_row`, `update_director_view`, `store_position_view` names used consistently across tasks.

### Task 2.5 (added after T2 benchmark): columnar ghost POSITION update + FORCE reduction

**Measured motivation:** post-T2 ratios lj4 1.055 / p3m1 1.056 / p3m4 1.070 — residue isolated to the per-step ghost path (per-particle value serialization from strided columns). Implements spec section 2's ghost design (per-column gathers; contiguous ghost row ranges) for the two per-step field groups.

**Files:** src/core/ghosts.cpp (+ CellStructure access to store views).

**Design:** when `data_parts == GHOSTTRANS_POSITION` (the per-step ghost update) or `== GHOSTTRANS_FORCE` (reduction), take a columnar fast path: for each part_list, rows form a contiguous range `[first_row, first_row + size)` (cell-sorted store; ASSERT contiguity in debug); pack/unpack loops read/write the position/image/quaternion (or force/torque) columns directly by row — no Particle accessors, no proxies. Wire layout stays byte-identical (per-particle field order, same sizes). SAVE+shift path applies `+shift` and `fold_position` inside the bulk loop. Any other data_parts combination falls back to the existing per-particle `serialize_and_reduce` path (resort-time comms are not hot). GHOST_LOCL cell_cell_transfer gets the same bulk treatment. Identity must stay bitwise; the full ghost-mode battery (send/recv, BCST, RDCE, LOCL; 1/4 ranks; lees_edwards shift; hybrid) gates it.
