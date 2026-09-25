# ParticleStore Migration — Phase 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `ParticleStore` real row bookkeeping (cell-sorted rebuild, dirty-marking, `store_row` in `Particle`) and evict force + torque out of the `Particle` struct into Kokkos columns as the single source of truth.

**Architecture:** Force/torque become `LayoutLeft` Kokkos columns in `ParticleStore` (owned by `CellStructure`). Every stored particle carries a transitional back-pointer + row index; a lazy, purely rank-local `synchronize` rebuilds rows in cell-traversal order whenever topology changed. `Particle::force()` returns a write-through `VectorReference` proxy (non-const) or a `Utils::Vector3d` value (const). The struct members are deleted in the same commit that flips the source of truth (spec section 4, single ownership). Ghost force reduction keeps its wire format but serializes force as explicit values, never via the struct. The Python `f`/`torque_lab` getters switch from the fetch-cache detached copy to an owner-side MPI fetch.

**Tech Stack:** C++20, Kokkos (Views, DualView, ScatterView already in place), Boost.MPI/Boost.Test, Python test suite, benchmark gate from phase 0.

**Spec:** `docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md` (section 3 phase 2, sections 4, 5).

## Global Constraints

- Worktree `/tikhome/weeber/es/.claude/worktrees/eliminate_Particle` (branch `worktree-eliminate_Particle`) only; all commands from the worktree root. Never touch other worktrees.
- `make -j8` only; never `-j$(nproc)`. Do not re-run cmake configure (build/ exists).
- Tests green at every commit. Phase gate: full Python suite green; benchmarks (LJ, P3M at `--particles_per_core 1000` on 1 and 4 MPI ranks; `--particles_per_core 4000` with 4 OMP threads) at most 5% cumulative regression vs `maintainer/benchmarks/baselines/phase0-baseline.csv`, measured with `maintainer/benchmarks/benchmark_gate.py` (busy check = foreign CPU; never bypass on a loaded machine — the machine is shared).
- **Numerical identity** (spec section 5): a fixed-seed serial LJ and P3M trajectory must be bitwise-identical in positions and forces before and after this phase (Task 1 records the reference; later tasks compare).
- **Single ownership** (spec section 4): the flip commit removes `ParticleForce` members from `Particle` in the same commit that makes columns authoritative. No transitional dual storage.
- Pure Kokkos storage; `LayoutLeft` (component-major) for vector columns per spec.
- Full-word identifier naming.
- CUDA stays OFF (CPU-only build); code under `ESPRESSO_CUDA`/GPU ifdefs must still be updated consistently but cannot be compile-tested here — keep those edits mechanical (proxy `operator[]` compatible).
- Pre-commit hook runs clang-format; if it reformats, re-stage and re-commit. Commit messages end with:
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- Known machinery (do not redesign): `m_local_force`/`m_scatter_force` scratch accumulation and the contribute-back loop in `src/core/forces.cpp:209-262` stay as they are; after the flip they write through the new accessors.

---

### Task 1: Numerical-identity reference capture

**Files:**
- Create: `maintainer/benchmarks/trajectory_identity.py`
- Create (uncommitted session artifact): `.superpowers/sdd/phase1-identity-reference.txt`

**Interfaces:**
- Produces: `./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj|p3m` printing exactly one line `<mode> <sha256-of-positions-and-forces>`; the reference file with both lines. Tasks 4, 5, 7, 8 re-run and diff against it.

- [ ] **Step 1: Write the script**

```python
#!/usr/bin/env python3
#
# (GPL header block: copy from maintainer/benchmarks/benchmark_gate.py lines 1-18)
#
"""
Fixed-seed, serial, deterministic short trajectory for the ParticleStore
migration's bitwise-identity gate (spec section 5, phases 2-6).

Prints one line: '<mode> <sha256>' where the hash covers the bitwise
content of positions and forces after the run. Any change in this hash
between two commits of the same phase means the storage migration
altered numerics and is a bug (phases 2-6 are pure storage moves).
"""

import argparse
import hashlib

import numpy as np

import espressomd
import espressomd.electrostatics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mode", choices=["lj", "p3m"], required=True)
args = parser.parse_args()

np.random.seed(42)

system = espressomd.System(box_l=[12.0, 12.0, 12.0])
system.time_step = 0.001
system.cell_system.skin = 0.4

# 6x6x6 lattice with a small deterministic perturbation
positions = []
for i in range(6):
    for j in range(6):
        for k in range(6):
            positions.append([2.0 * i + 1.0, 2.0 * j + 1.0, 2.0 * k + 1.0])
positions = np.array(positions) + 0.1 * np.random.random((216, 3))

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0, cutoff=2.5, shift="auto")

if args.mode == "lj":
    system.part.add(pos=positions)
else:
    charges = np.resize([1.0, -1.0], 216)
    system.part.add(pos=positions, q=charges)
    # fully pinned parameters: no tuning, fully deterministic
    solver = espressomd.electrostatics.P3M(
        prefactor=1.0, accuracy=1e-4, mesh=[16, 16, 16], cao=6,
        r_cut=3.5, alpha=0.85, tune=False)
    system.electrostatics.solver = solver

system.integrator.run(20)

digest = hashlib.sha256()
digest.update(np.copy(system.part.all().pos).tobytes())
digest.update(np.copy(system.part.all().f).tobytes())
print(f"{args.mode} {digest.hexdigest()}")
```

- [ ] **Step 2: Run both modes twice — verify determinism**

```bash
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode p3m
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode p3m
```
Expected: each mode prints the SAME hash on both runs. If not deterministic (e.g. P3M node-grid or FFT nondeterminism), investigate before proceeding — the gate is worthless without run-to-run determinism. If P3M cannot be made deterministic with pinned parameters, report BLOCKED.

- [ ] **Step 3: Record the reference**

```bash
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj > .superpowers/sdd/phase1-identity-reference.txt
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode p3m >> .superpowers/sdd/phase1-identity-reference.txt
cat .superpowers/sdd/phase1-identity-reference.txt
```
Expected: two lines. This file is the phase-1 reference (valid for this machine + build flags only; intentionally not committed).

- [ ] **Step 4: Commit the script only**

```bash
git add maintainer/benchmarks/trajectory_identity.py
git commit -m "benchmarks: add fixed-seed trajectory identity gate (spec section 5)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: VectorReference proxy

**Files:**
- Create: `src/core/particle_store/VectorReference.hpp`
- Create: `src/core/unit_tests/VectorReference_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Produces (used verbatim by Tasks 3 and 7): `class VectorReference` — write-through view of one particle's 3-vector inside a component-major column. Constructed from `(double *base, std::size_t stride)` where component `j` lives at `base[j * stride]`. Operations: implicit `operator Utils::Vector3d() const`, `operator=(Utils::Vector3d const &)`, `operator+=`, `operator-=`, `operator*=(double)`, `operator[](std::size_t) -> double &` (+ const overload), `norm()`, `norm2()`.

- [ ] **Step 1: Write the failing test**

`src/core/unit_tests/VectorReference_test.cpp` (GPL header from `Particle_test.cpp` lines 1-18, then):

```cpp
#define BOOST_TEST_MODULE VectorReference test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>

// simulate a component-major (LayoutLeft) column with 4 rows:
// memory is [x0 x1 x2 x3 | y0 y1 y2 y3 | z0 z1 z2 z3], stride 4
struct ColumnFixture {
  std::array<double, 12> storage{};
  static constexpr std::size_t stride = 4;
  VectorReference row(std::size_t i) {
    return VectorReference(storage.data() + i, stride);
  }
};

BOOST_FIXTURE_TEST_CASE(assignment_writes_through, ColumnFixture) {
  row(1) = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(storage[1], 1.);  // x1
  BOOST_CHECK_EQUAL(storage[5], 2.);  // y1
  BOOST_CHECK_EQUAL(storage[9], 3.);  // z1
}

BOOST_FIXTURE_TEST_CASE(conversion_reads_components, ColumnFixture) {
  storage[2] = 4.;  // x2
  storage[6] = 5.;  // y2
  storage[10] = 6.; // z2
  Utils::Vector3d const value = row(2);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 5.);
  BOOST_CHECK_EQUAL(value[2], 6.);
}

BOOST_FIXTURE_TEST_CASE(compound_operators, ColumnFixture) {
  row(0) = Utils::Vector3d{1., 1., 1.};
  row(0) += Utils::Vector3d{1., 2., 3.};
  row(0) -= Utils::Vector3d{0., 1., 0.};
  row(0) *= 2.;
  Utils::Vector3d const value = row(0);
  BOOST_CHECK_EQUAL(value[0], 4.);
  BOOST_CHECK_EQUAL(value[1], 4.);
  BOOST_CHECK_EQUAL(value[2], 8.);
}

BOOST_FIXTURE_TEST_CASE(subscript_and_norms, ColumnFixture) {
  row(3)[0] = 3.;
  row(3)[1] = 0.;
  row(3)[2] = 4.;
  BOOST_CHECK_EQUAL(row(3)[2], 4.);
  BOOST_CHECK_EQUAL(row(3).norm2(), 25.);
  BOOST_CHECK_EQUAL(row(3).norm(), 5.);
}

BOOST_FIXTURE_TEST_CASE(proxy_copies_alias_the_same_row, ColumnFixture) {
  auto reference = row(1);
  auto alias = reference; // copying the proxy copies the reference, not data
  alias = Utils::Vector3d{7., 8., 9.};
  BOOST_CHECK_EQUAL(storage[1], 7.);
}
```

Register in `src/core/unit_tests/CMakeLists.txt` next to the `ParticleStore_test.cpp` entry:

```cmake
espresso_unit_test(SRC VectorReference_test.cpp DEPENDS espresso::utils)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `make -C build -j8 VectorReference_test 2>&1 | tail -5`
Expected: FAIL — `particle_store/VectorReference.hpp: No such file or directory`

- [ ] **Step 3: Write the header**

`src/core/particle_store/VectorReference.hpp` (GPL header, then):

```cpp
#pragma once

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>

/**
 * @brief Write-through reference to one particle's 3-vector stored in a
 * component-major (LayoutLeft) column.
 *
 * Component @c j of the referenced vector lives at <tt>base[j * stride]</tt>.
 * The proxy is a cheap value type (pointer + stride); copying it copies the
 * reference, not the data. Reads convert to @ref Utils::Vector3d by value.
 *
 * Note: because @ref Utils::Vector arithmetic operators are function
 * templates, a VectorReference does not participate in template argument
 * deduction — expressions like <tt>proxy - vector</tt> require an explicit
 * conversion: <tt>Utils::Vector3d(proxy) - vector</tt>.
 */
class VectorReference {
  double *m_base;
  std::size_t m_stride;

public:
  VectorReference(double *base, std::size_t stride)
      : m_base(base), m_stride(stride) {
    assert(base != nullptr);
  }

  operator Utils::Vector3d() const {
    return {m_base[0u], m_base[m_stride], m_base[2u * m_stride]};
  }

  VectorReference &operator=(Utils::Vector3d const &value) {
    m_base[0u] = value[0u];
    m_base[m_stride] = value[1u];
    m_base[2u * m_stride] = value[2u];
    return *this;
  }

  VectorReference &operator+=(Utils::Vector3d const &value) {
    m_base[0u] += value[0u];
    m_base[m_stride] += value[1u];
    m_base[2u * m_stride] += value[2u];
    return *this;
  }

  VectorReference &operator-=(Utils::Vector3d const &value) {
    m_base[0u] -= value[0u];
    m_base[m_stride] -= value[1u];
    m_base[2u * m_stride] -= value[2u];
    return *this;
  }

  VectorReference &operator*=(double const factor) {
    m_base[0u] *= factor;
    m_base[m_stride] *= factor;
    m_base[2u * m_stride] *= factor;
    return *this;
  }

  double &operator[](std::size_t const j) {
    assert(j < 3u);
    return m_base[j * m_stride];
  }
  double const &operator[](std::size_t const j) const {
    assert(j < 3u);
    return m_base[j * m_stride];
  }

  double norm2() const {
    auto const x = m_base[0u];
    auto const y = m_base[m_stride];
    auto const z = m_base[2u * m_stride];
    return x * x + y * y + z * z;
  }
  double norm() const { return std::sqrt(norm2()); }
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `make -C build -j8 VectorReference_test && ctest --test-dir build -R VectorReference_test --output-on-failure`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/core/particle_store/VectorReference.hpp src/core/unit_tests/VectorReference_test.cpp src/core/unit_tests/CMakeLists.txt
git commit -m "core: add VectorReference write-through proxy (migration phase 2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: ParticleStore columns, rebuild machinery, CellStructure integration

**Files:**
- Modify: `src/core/particle_store/ParticleStore.hpp` (grow the phase-0 skeleton)
- Modify: `src/core/Particle.hpp` (transitional attach members — accessors NOT yet flipped)
- Modify: `src/core/cell_system/CellStructure.hpp` + `src/core/cell_system/CellStructure.cpp` (member, dirty marks, synchronize)
- Modify: `src/core/unit_tests/ParticleStore_test.cpp` (extend)

**Interfaces:**
- Consumes: `VectorReference(double*, std::size_t)` (Task 2).
- Produces (used verbatim by Tasks 6-7):
  - `ParticleStore::mark_dirty()`, `bool is_dirty() const`
  - `ParticleStore::begin_rebuild(std::size_t n_local, std::size_t n_ghost)` — snapshots old columns, allocates fresh zero-initialized columns sized `n_local + n_ghost`
  - `ParticleStore::assign_row(Particle &p, int row)` — copies the particle's old force/torque values into the new row when it had a valid old row, then attaches the particle (`p.attach_to_store(*this, row)`)
  - `ParticleStore::finish_rebuild()` — drops old columns, clears dirty
  - `VectorReference ParticleStore::force_reference(int row)`, `Utils::Vector3d force_value(int row) const`; torque analogues under `ESPRESSO_ROTATION`
  - `Particle::attach_to_store(ParticleStore &, int)`, `Particle::store_row() const -> int` (`-1` = detached), `ParticleStore *Particle::store() const`
  - `CellStructure::particle_store() -> ParticleStore &`, `CellStructure::mark_particle_store_dirty()`, `CellStructure::ensure_particle_store_synchronized()` (O(1) when clean; rebuild is purely rank-local — safe to call from any rank independently, no MPI)

- [ ] **Step 1: Write the failing tests** (append to `src/core/unit_tests/ParticleStore_test.cpp`; add includes `#include "Particle.hpp"` and `#include <utils/Vector.hpp>`)

```cpp
BOOST_AUTO_TEST_CASE(rebuild_assigns_rows_and_zero_initializes) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.mark_dirty();
  BOOST_CHECK(store.is_dirty());
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  BOOST_CHECK(not store.is_dirty());
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 2ul);
  BOOST_CHECK_EQUAL(p0.store_row(), 0);
  BOOST_CHECK_EQUAL(p1.store_row(), 1);
  Utils::Vector3d const f0 = store.force_value(0);
  BOOST_CHECK_EQUAL(f0.norm2(), 0.);
}

BOOST_AUTO_TEST_CASE(rebuild_preserves_values_by_old_row) {
  ParticleStore store{};
  Particle p0{}, p1{};
  store.begin_rebuild(2u, 0u);
  store.assign_row(p0, 0);
  store.assign_row(p1, 1);
  store.finish_rebuild();
  store.force_reference(p0.store_row()) = Utils::Vector3d{1., 2., 3.};
  store.force_reference(p1.store_row()) = Utils::Vector3d{4., 5., 6.};

  // simulate a resort that swaps the two particles' order and adds one
  Particle p2{};
  store.mark_dirty();
  store.begin_rebuild(3u, 0u);
  store.assign_row(p1, 0);
  store.assign_row(p0, 1);
  store.assign_row(p2, 2);
  store.finish_rebuild();

  Utils::Vector3d const f_p1 = store.force_value(p1.store_row());
  Utils::Vector3d const f_p0 = store.force_value(p0.store_row());
  Utils::Vector3d const f_p2 = store.force_value(p2.store_row());
  BOOST_CHECK_EQUAL(f_p1[0], 4.); // p1 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p0[2], 3.); // p0 kept its values at its new row
  BOOST_CHECK_EQUAL(f_p2.norm2(), 0.); // new particle zero-initialized
}

BOOST_AUTO_TEST_CASE(ghost_rows_follow_locals) {
  ParticleStore store{};
  Particle p_local{}, p_ghost{};
  store.begin_rebuild(1u, 1u);
  store.assign_row(p_local, 0);
  store.assign_row(p_ghost, 1);
  store.finish_rebuild();
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 1ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 1ul);
}
```

- [ ] **Step 2: Run to verify failure**

Run: `make -C build -j8 ParticleStore_test 2>&1 | tail -5`
Expected: FAIL (missing `mark_dirty` etc.)

- [ ] **Step 3: Implement `ParticleStore`**

Replace the body of `src/core/particle_store/ParticleStore.hpp` (keep GPL header and the class doc comment, update the "scaffolding" note):

```cpp
#pragma once

#include <config/config.hpp>

#include "particle_store/VectorReference.hpp"

#include <utils/Vector.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

#include <cassert>
#include <cstddef>

class Particle; // attach_to_store is defined in Particle.hpp

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns per-particle quantities in a single index space: local particles
 * first (cell-traversal order), ghosts appended. Fields are component-major
 * (LayoutLeft) Kokkos columns; see
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 *
 * Migration phase 2: force and torque columns (observables). Rebuild
 * protocol: mark_dirty() on any topology change; the owner (CellStructure)
 * later runs begin_rebuild / assign_row-per-particle / finish_rebuild.
 * Rebuild preserves values by old row and zero-initializes new rows.
 * Rebuilds are purely rank-local (no MPI).
 */
class ParticleStore {
public:
  using Column = Kokkos::DualView<double *[3], Kokkos::LayoutLeft>;

  std::size_t number_of_local_particles() const {
    return m_number_of_local_particles;
  }
  std::size_t number_of_ghost_particles() const {
    return m_number_of_ghost_particles;
  }
  std::size_t number_of_particles() const {
    return m_number_of_local_particles + m_number_of_ghost_particles;
  }

  void mark_dirty() { m_dirty = true; }
  bool is_dirty() const { return m_dirty; }

  void begin_rebuild(std::size_t number_of_local_particles,
                     std::size_t number_of_ghost_particles);
  void assign_row(Particle &particle, int row);
  void finish_rebuild();

  VectorReference force_reference(int const row) {
    return column_reference(m_force, row);
  }
  Utils::Vector3d force_value(int const row) const {
    return column_value(m_force, row);
  }
#ifdef ESPRESSO_ROTATION
  VectorReference torque_reference(int const row) {
    return column_reference(m_torque, row);
  }
  Utils::Vector3d torque_value(int const row) const {
    return column_value(m_torque, row);
  }
#endif

private:
  VectorReference column_reference(Column &column, int const row) {
    assert(row >= 0 and
           static_cast<std::size_t>(row) < number_of_particles());
    auto &view = column.view_host();
    return VectorReference(view.data() + row, view.stride_1());
  }
  Utils::Vector3d column_value(Column const &column, int const row) const {
    assert(row >= 0 and
           static_cast<std::size_t>(row) < number_of_particles());
    auto const &view = column.view_host();
    return {view(row, 0), view(row, 1), view(row, 2)};
  }

  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
  bool m_dirty = false;
  Column m_force;
#ifdef ESPRESSO_ROTATION
  Column m_torque;
#endif
  // previous-generation columns, alive between begin_rebuild/finish_rebuild
  Column m_old_force;
#ifdef ESPRESSO_ROTATION
  Column m_old_torque;
#endif
  std::size_t m_old_number_of_particles = 0u;
};
```

Implementation of the three rebuild functions — because `assign_row` needs the complete `Particle` type, put them in a new small implementation file `src/core/particle_store/ParticleStore.cpp` (add it to the core sources: find where `cells.cpp` is listed in `src/core/CMakeLists.txt` and add `particle_store/ParticleStore.cpp` alongside):

```cpp
// GPL header
#include "particle_store/ParticleStore.hpp"

#include "Particle.hpp"

void ParticleStore::begin_rebuild(
    std::size_t const number_of_local_particles,
    std::size_t const number_of_ghost_particles) {
  m_old_force = m_force;
#ifdef ESPRESSO_ROTATION
  m_old_torque = m_torque;
#endif
  m_old_number_of_particles =
      m_number_of_local_particles + m_number_of_ghost_particles;
  m_number_of_local_particles = number_of_local_particles;
  m_number_of_ghost_particles = number_of_ghost_particles;
  auto const total = number_of_local_particles + number_of_ghost_particles;
  // fresh zero-initialized allocation (Kokkos default-initializes to zero)
  m_force = Column("particle_store::force", total);
#ifdef ESPRESSO_ROTATION
  m_torque = Column("particle_store::torque", total);
#endif
}

void ParticleStore::assign_row(Particle &particle, int const row) {
  auto const old_row = particle.store_row();
  if (particle.store() == this and old_row >= 0 and
      static_cast<std::size_t>(old_row) < m_old_number_of_particles) {
    auto const &old_force = m_old_force.view_host();
    auto &new_force = m_force.view_host();
    for (std::size_t j = 0u; j < 3u; ++j) {
      new_force(row, j) = old_force(old_row, j);
    }
#ifdef ESPRESSO_ROTATION
    auto const &old_torque = m_old_torque.view_host();
    auto &new_torque = m_torque.view_host();
    for (std::size_t j = 0u; j < 3u; ++j) {
      new_torque(row, j) = old_torque(old_row, j);
    }
#endif
  }
  particle.attach_to_store(*this, row);
}

void ParticleStore::finish_rebuild() {
  m_old_force = Column{};
#ifdef ESPRESSO_ROTATION
  m_old_torque = Column{};
#endif
  m_old_number_of_particles = 0u;
  m_dirty = false;
}
```

- [ ] **Step 4: Add transitional attach members to `Particle`**

In `src/core/Particle.hpp`, forward-declare `class ParticleStore;` before `struct Particle`, and inside `struct Particle`'s private section add (next to the other members):

```cpp
  /** Transitional (migration phase 2): row of this particle in the
   *  ParticleStore, -1 while detached. Rank-local; never serialized. */
  ParticleStore *m_particle_store = nullptr;
  int m_store_row = -1;
```

and in the public section:

```cpp
  void attach_to_store(ParticleStore &store, int const row) {
    m_particle_store = &store;
    m_store_row = row;
  }
  auto store() const { return m_particle_store; }
  auto store_row() const { return m_store_row; }
```

Do NOT touch `Particle::serialize` (the members stay rank-local: deserialized particles arrive detached, which is correct).

- [ ] **Step 5: CellStructure integration**

In `src/core/cell_system/CellStructure.hpp`: `#include "particle_store/ParticleStore.hpp"`, add member `ParticleStore m_particle_store;` and public API:

```cpp
  auto &particle_store() { return m_particle_store; }
  void mark_particle_store_dirty() { m_particle_store.mark_dirty(); }
  /** @brief Rebuild the store row assignment if topology changed.
   *  Purely rank-local; O(1) when the store is clean. */
  void ensure_particle_store_synchronized();
```

In `src/core/cell_system/CellStructure.cpp` implement:

```cpp
void CellStructure::ensure_particle_store_synchronized() {
  if (not m_particle_store.is_dirty()) {
    return;
  }
  auto const n_local = count_local_particles();
  std::size_t n_ghost = 0u;
  for (auto const &p : ghost_particles()) {
    static_cast<void>(p);
    ++n_ghost;
  }
  m_particle_store.begin_rebuild(n_local, n_ghost);
  int row = 0;
  for (auto &p : local_particles()) {
    m_particle_store.assign_row(p, row++);
  }
  for (auto &p : ghost_particles()) {
    m_particle_store.assign_row(p, row++);
  }
  m_particle_store.finish_rebuild();
}
```

(Check how `count_local_particles()` is spelled in this class — use the existing helper; if none fits, count `local_particles()` the same way as ghosts.)

Add `mark_particle_store_dirty();` to these methods (all in CellStructure.cpp/.hpp): `add_particle`, `add_local_particle`, `remove_particle`, `remove_all_particles`, `resort_particles` (after `m_decomposition->resort(...)` returns), `ghosts_count` (after the `GHOSTTRANS_PARTNUM` communication), and `set_particle_decomposition`. Also mark it initially dirty on construction so the first synchronize builds rows (the member initializer `bool m_dirty = false;` in ParticleStore stays; instead call `mark_particle_store_dirty()` at the end of the CellStructure constructor(s)).

- [ ] **Step 6: Build, run tests**

Run: `make -C build -j8 && ctest --test-dir build -R "ParticleStore_test|VectorReference_test" --output-on-failure && ctest --test-dir build -L unit_test --output-on-failure 2>&1 | tail -3`
Expected: all pass. (The store is wired but nothing reads it yet — force still lives in the struct. That is intentional staging: Task 7 flips the source of truth. Note for the reviewer: the dead-until-Task-7 machinery is the point of this task's gate — it proves rebuild bookkeeping independently before the risky flip.)

- [ ] **Step 7: Quick Python sanity + commit**

Run: `ctest --test-dir build -R "^particle$|cellsystem" --output-on-failure`
Expected: pass.

```bash
git add src/core/particle_store/ src/core/Particle.hpp src/core/cell_system/CellStructure.hpp src/core/cell_system/CellStructure.cpp src/core/CMakeLists.txt src/core/unit_tests/ParticleStore_test.cpp
git commit -m "core: ParticleStore row bookkeeping and force/torque columns (phase 2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Call-site preparation audit (no behavior change)

**Files:**
- Modify: `src/core/electrostatics/elc.cpp:617` and `:793` (reference bindings)
- Modify: `src/core/magnetostatics/dp3m_heffte.impl.hpp:896` (reference binding)
- Modify: `src/core/forces.cpp` (`init_forces_and_thermostat`, `ghosts_reset_forces` call target — see below)
- Modify: `src/core/stokesian_dynamics/sd_interface.cpp:52-53`
- Modify: `src/core/unit_tests/rotation_test.cpp:179`
- Possibly more sites found by the greps below.

**Interfaces:**
- Consumes: nothing new — every edit keeps today's `Utils::Vector3d&` accessors working identically.
- Produces: a codebase where `p.force()`/`p.torque()` appear ONLY in proxy-compatible patterns: `= v`, `+= v`, `-= v`, `*= s`, `[j]`, `.norm()`, `.norm2()`, passed as `Utils::Vector3d const&`, or explicitly converted. All uses of `force_and_torque()` outside ghosts.cpp are eliminated (ghosts.cpp is Task 5's job).

- [ ] **Step 1: Fix the three reference bindings**

In `src/core/electrostatics/elc.cpp` (~617 and ~793), pattern:

```cpp
    auto &force = p.force();
    ... force[i] += ...;
```

becomes:

```cpp
    Utils::Vector3d force = p.force();
    ... force[i] += ...;   // (loop body unchanged)
    p.force() = force;     // write back after the last component update
```

Apply the same shape at `src/core/magnetostatics/dp3m_heffte.impl.hpp:896` for `torque`.

- [ ] **Step 2: Eliminate `force_and_torque()` uses outside ghosts.cpp**

Greps to enumerate (verify — line numbers may have drifted):

```bash
grep -rn "force_and_torque" src/core src/script_interface --include=*.cpp --include=*.hpp | grep -v ghosts.cpp
```

Expected sites and fixes:
- `src/core/forces.cpp:113` — `p.force_and_torque() = external_force(p);` → 
  ```cpp
  auto const external = external_force(p);
  p.force() = external.f;
  #ifdef ESPRESSO_ROTATION
  p.torque() = external.torque;
  #endif
  ```
  (`external_force` returns `ParticleForce`; check its actual name/return at the call site and keep the `ROTATION` guard consistent with `ParticleForce`'s definition.)
- `CellStructure::ghosts_reset_forces` (find via `grep -rn "ghosts_reset_forces" src/core`) — if it does `p.force_and_torque() = {}` (or `= ParticleForce{}`), replace with `p.force() = {};` plus `p.torque() = {};` under `ESPRESSO_ROTATION`.
- `src/core/stokesian_dynamics/sd_interface.cpp:52-53` — constructor member init `ext_force(p.force_and_torque())` → `ext_force(ParticleForce{p.force()
#ifdef ESPRESSO_ROTATION
      , p.torque()
#endif
  })` — check `ParticleForce`'s constructors (it has `(f)` and `(f, torque)`); pick the matching one.
- `Particle::force_and_torque()` accessor itself stays for now (deleted in Task 7).

- [ ] **Step 3: Fix proxy-hostile read patterns**

```bash
grep -rn "auto [a-z_]* = p[a-z_.]*\.force()\|auto [a-z_]* = p[a-z_.]*\.torque()\|auto const [a-z_]* = [a-z_.]*\.torque()\|auto const [a-z_]* = [a-z_.]*\.force()" src/core src/script_interface
grep -rn "\.force() [-+*]\|\.torque() [-+*]" src/core src/script_interface testsuite 2>/dev/null | grep -v "+=\|-=\|*="
```

For every hit: `auto x = p.force();` → `Utils::Vector3d x = p.force();` and binary arithmetic `p.force() - v` → `Utils::Vector3d(p.force()) - v`. Known: `src/core/unit_tests/rotation_test.cpp:179` (`auto const t_out = p.torque();` → `Utils::Vector3d const t_out = p.torque();`) and `src/core/unit_tests/lb_particle_coupling_test.cpp:343` (`(p.force() - expected).norm()` → `(Utils::Vector3d(p.force()) - expected).norm()`). Fix all hits the same way.

- [ ] **Step 4: Build, test, identity check**

Run: `make -C build -j8 && ctest --test-dir build -L unit_test --output-on-failure 2>&1 | tail -3`
Expected: pass.
Run: `./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj && ./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode p3m` and compare with `.superpowers/sdd/phase1-identity-reference.txt`.
Expected: hashes IDENTICAL (these edits reorder no arithmetic — `force = p.force(); ...; p.force() = force` performs the same operations on the same values).
Run: `ctest --test-dir build -R "elc|dipol|magneto|stokesian|rotation|coulomb" --output-on-failure`
Expected: pass (covers the touched solvers).

- [ ] **Step 5: Commit**

```bash
git add -A src/core src/script_interface
git commit -m "core: prepare force/torque call sites for proxy accessors (phase 2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Ghost FORCE path — value-based field serialization

**Files:**
- Modify: `src/core/ghosts.cpp` (`serialize_and_reduce` FORCE branch, `calc_transmit_size`, `add_forces_from_recv_buffer`)

**Interfaces:**
- Consumes: today's struct-backed accessors (this task lands BEFORE the flip; behavior and wire format are unchanged).
- Produces: a ghost FORCE path that never serializes `Particle` structs or `ParticleForce` members through particle memory: SAVE writes explicit `Utils::Vector3d` values read via accessors; LOAD does `p.force() += value`; transmit size for FORCE is a compile-time constant; FORCE is asserted to never mix with other data parts in one communication.

- [ ] **Step 1: Rework `serialize_and_reduce`'s FORCE branch** (src/core/ghosts.cpp:~197-215)

Replace the current branch (which does `ar & p.force();` on SAVE) with value copies in BOTH directions:

```cpp
  if (data_parts & GHOSTTRANS_FORCE) {
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d force;
      ar & force;
      p.force() += force;
    } else {
      Utils::Vector3d force = p.force();
      ar & force;
    }
#ifdef ESPRESSO_ROTATION
    if (policy == ReductionPolicy::UPDATE and
        direction == SerializationDirection::LOAD) {
      Utils::Vector3d torque;
      ar & torque;
      p.torque() += torque;
    } else {
      Utils::Vector3d torque = p.torque();
      ar & torque;
    }
#endif
  }
```

(The LOAD+UPDATE side is already value-based today; the change is the SAVE side and the non-UPDATE LOAD side — grep the existing branch carefully and keep MOVE-policy LOAD semantics identical: if the current code does `ar & p.force()` for MOVE+LOAD, replace with `Utils::Vector3d force; ar & force; p.force() = force;`.)

- [ ] **Step 2: Constant-size FORCE in `calc_transmit_size`** (src/core/ghosts.cpp:~231-251)

At the top of `calc_transmit_size(BoxGeometry const&, unsigned data_parts)` add:

```cpp
  if (data_parts & GHOSTTRANS_FORCE) {
    // Forces are always communicated alone (collect_ghost_force_comm).
    assert(data_parts == GHOSTTRANS_FORCE);
#ifdef ESPRESSO_ROTATION
    return 6ul * sizeof(double);
#else
    return 3ul * sizeof(double);
#endif
  }
```

This removes the need to serialize a default-constructed `Particle{}` for the FORCE part (which becomes impossible after the flip — a detached particle has no force row).

- [ ] **Step 3: Rework `add_forces_from_recv_buffer`** (src/core/ghosts.cpp:~356-367)

Replace the `ParticleForce pf; archiver >> pf; part.force_and_torque() += pf;` body with:

```cpp
      Utils::Vector3d force;
      archiver >> force;
      part.force() += force;
#ifdef ESPRESSO_ROTATION
      Utils::Vector3d torque;
      archiver >> torque;
      part.torque() += torque;
#endif
```

(Wire layout is identical: `ParticleForce` serializes as f then torque.)

- [ ] **Step 4: Build, test (MPI-heavy), identity check**

Run: `make -C build -j8 && ctest --test-dir build -L unit_test --output-on-failure 2>&1 | tail -3 && ctest --test-dir build -R "particle|cellsystem|coulomb|lees_edwards" --output-on-failure`
Expected: pass (4-rank variants exercise GHOST_RDCE, GHOST_LOCL, and send/recv force reduction).
Identity: run both modes, compare to reference — IDENTICAL (same values, same order, same wire bytes).

- [ ] **Step 5: Commit**

```bash
git add src/core/ghosts.cpp
git commit -m "core: ghost force reduction via explicit field values (phase 2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Python force/torque getters — owner-side fetch

**Files:**
- Modify: `src/core/particle_node.hpp`, `src/core/particle_node.cpp`
- Modify: `src/script_interface/particle_data/ParticleHandle.cpp:253` (f getter) and `:402-405` (torque_lab getter)

**Interfaces:**
- Consumes: nothing from the store yet (pre-flip: owner reads struct-backed accessors on the LIVE particle).
- Produces: `Utils::Vector3d get_particle_force(int p_id)` and (under `ESPRESSO_ROTATION`) `Utils::Vector3d get_particle_torque_lab(int p_id)` in particle_node — MPI collective fetches that read from the owning rank's live particle, never from the fetch-cache detached copy. After the flip these are the ONLY remote read paths for force/torque.

- [ ] **Step 1: Add the fetch functions**

`src/core/particle_node.hpp` — declare:

```cpp
/** @brief Get the force on a particle, fetched from the owning rank. */
Utils::Vector3d get_particle_force(int p_id);
#ifdef ESPRESSO_ROTATION
/** @brief Get a particle's torque in the lab frame, from the owning rank. */
Utils::Vector3d get_particle_torque_lab(int p_id);
#endif
```

`src/core/particle_node.cpp` — implement following the existing callback pattern of `mpi_send_particle_data_local` (line ~159; note `REGISTER_CALLBACK` and `Communication::mpiCallbacks().call_all(...)` usage — copy its structure exactly, including the found-count assertion):

```cpp
static auto get_local_particle_property(
    int p_id, Utils::Vector3d (*getter)(Particle const &)) {
  auto const p = get_cell_structure().get_local_particle(p_id);
  auto const found = (p != nullptr) and not p->is_ghost();
  assert(1 == boost::mpi::all_reduce(::comm_cart, static_cast<int>(found),
                                     std::plus<>()) &&
         "particle not found exactly once");
  auto const local_value = found ? getter(*p) : Utils::Vector3d{};
  return boost::mpi::all_reduce(::comm_cart, local_value, std::plus<>());
}

static void mpi_get_particle_force_local(int p_id) {
  get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}

REGISTER_CALLBACK(mpi_get_particle_force_local)

Utils::Vector3d get_particle_force(int p_id) {
  Communication::mpiCallbacks().call_all(mpi_get_particle_force_local, p_id);
  return get_local_particle_property(
      p_id, [](Particle const &p) { return Utils::Vector3d(p.force()); });
}
```

CAUTION: `call_all` runs the callback on the worker ranks while the head rank executes the follow-up call itself — mirror exactly how `get_particle_data` at line ~172-192 pairs `call_all` with head-node work, so every rank participates in both `all_reduce`s exactly once. If `call_all` also runs on the head node in this codebase's `MpiCallbacks` semantics (check `get_particle_data`!), adjust so the reduction executes once per rank — the existing function is the authoritative template. Non-capturing lambdas convert to the function pointer.

`get_particle_torque_lab` is identical except the getter does the frame conversion ON the owner (include `"rotation.hpp"`):

```cpp
      [](Particle const &p) {
        return convert_vector_body_to_space(p, Utils::Vector3d(p.torque()));
      }
```

- [ ] **Step 2: Point ParticleHandle at the new fetch**

`src/script_interface/particle_data/ParticleHandle.cpp`:
- line ~253: `[this]() { return get_particle_data(m_pid).force(); }` → `[this]() { return get_particle_force(m_pid); }`
- torque_lab getter (~402): 
  ```cpp
       [this]() {
         auto &p = get_particle_data(m_pid);
         return convert_vector_body_to_space(p, p.torque());
       }},
  ```
  → `[this]() { return get_particle_torque_lab(m_pid); }},`

- [ ] **Step 3: Build + targeted Python tests (multi-rank matters: remote getter path)**

Run: `make -C build -j8 && ctest --test-dir build -R "^particle$|particle_slice|rotation|langevin" --output-on-failure`
Expected: pass, including 4-rank variants (a particle owned by rank != 0 exercises the new fetch).

- [ ] **Step 4: Commit**

```bash
git add src/core/particle_node.hpp src/core/particle_node.cpp src/script_interface/particle_data/ParticleHandle.cpp
git commit -m "python: fetch particle force/torque from owning rank (phase 2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: THE FLIP — force/torque move to ParticleStore columns

**Files:**
- Modify: `src/core/Particle.hpp` (delete members, flip accessors, drop from serialize)
- Modify: `src/core/cell_system/CellStructure.cpp`/`integrate.cpp`/`forces.cpp` (sync points)
- Modify: `src/script_interface/particle_data/ParticleHandle.cpp` (`get_real_particle` sync)
- Modify: `src/core/particle_node.cpp` (sync in fetch callbacks)
- Modify: `src/core/unit_tests/rotation_test.cpp`, `src/core/unit_tests/lb_particle_coupling_test.cpp`, any other unit test that touches force/torque on hand-made particles (fixture below)
- Create: `src/core/unit_tests/ParticleStoreTestFixture.hpp`

**Interfaces:**
- Consumes: everything from Tasks 2-6.
- Produces: `Particle::force()` (non-const) returns `VectorReference`; `Particle::force() const` returns `Utils::Vector3d` BY VALUE; same for `torque()` under `ESPRESSO_ROTATION`; `force_and_torque()` is DELETED; `ParticleForce` remains as a standalone POD (still used by `external_force`, stokesian dynamics, and kernel plumbing) but is no longer a member of `Particle`.

- [ ] **Step 1: Flip the accessors and delete the members**

In `src/core/Particle.hpp`:
1. Delete the member `ParticleForce f;` from `struct Particle` and remove `ar & f;` from `Particle::serialize`. KEEP the `ParticleForce` struct definition (it is still a useful POD for `external_force()` returns and comm buffers).
2. Delete the `force_and_torque()` accessors (both overloads). Any compile error naming them is a Task-4 escape — fix that call site the Task-4 way.
3. Replace the accessors (add `#include "particle_store/ParticleStore.hpp"` — note this pulls Kokkos headers into Particle.hpp; that is accepted for the transition):

```cpp
  auto force() {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_reference(m_store_row);
  }
  Utils::Vector3d force() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->force_value(m_store_row);
  }
#ifdef ESPRESSO_ROTATION
  auto torque() {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_reference(m_store_row);
  }
  Utils::Vector3d torque() const {
    assert(m_particle_store != nullptr);
    return m_particle_store->torque_value(m_store_row);
  }
#endif
```

The const overloads returning values (not proxies) is deliberate: reads in const contexts get a real `Utils::Vector3d`, so template arithmetic keeps working.

- [ ] **Step 2: Wire the sync points** (each is one line: `ensure_particle_store_synchronized()` — O(1) when clean, rank-local, safe everywhere)

1. `src/core/integrate.cpp` — at the TOP of each integrator loop iteration (before `integrator_step_1`, which reads previous-step forces; mid-step particle creation by collision handling would otherwise leave new particles rowless when step 2 reads forces): `system.cell_structure->ensure_particle_store_synchronized();` — find the loop in `integrate(...)` (~line 585-706) and place it as the first statement of the loop body. ALSO place one immediately after the collision-handling call inside the loop if forces are read again before the loop iterates (check the order: if `handle_collisions` runs before `integrator_step_2`, the second sync is REQUIRED).
2. `src/core/forces.cpp` — first statement of `calculate_forces(...)`: same call (covers energies/pressures too if they route through it — they don't read forces, so forces.cpp is sufficient).
3. `src/script_interface/particle_data/ParticleHandle.cpp` — in `get_real_particle` (used by all setters): call `cell_structure.ensure_particle_store_synchronized();` before the lookup (it has a cell_structure reference — check the actual helper at ~lines 53-68 for how it reaches it).
4. `src/core/particle_node.cpp` — first statement of `get_local_particle_property` (Task 6): `get_cell_structure().ensure_particle_store_synchronized();`
5. `src/core/system/System.cpp` — find `on_observable_calc` (grep `on_observable_calc`); add the sync call so observables reading forces (TotalForce, ParticleForces) see valid rows. If no such hook exists, add the call at the top of `PidObservable::evaluate` in `src/core/observables/PidObservable.cpp` instead.

- [ ] **Step 3: Unit-test fixture for detached particles**

`src/core/unit_tests/ParticleStoreTestFixture.hpp` (GPL header, then):

```cpp
#pragma once

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"

#include <cstddef>

/** Attach hand-made particles to a standalone store so force/torque
 *  accessors work in unit tests (migration phase 2+). */
struct ParticleStoreTestFixture {
  ParticleStore store{};
  int next_row = 0;

  explicit ParticleStoreTestFixture(std::size_t capacity = 8u) {
    store.begin_rebuild(capacity, 0u);
    store.finish_rebuild();
  }
  void attach(Particle &p) { store.assign_row(p, next_row++); }
};
```

(`assign_row` outside begin/finish is fine here: it only copies old values when the particle was previously attached to the same store.) Adapt every unit test that compiles-or-asserts-fails: run `make -C build -j8 unit_tests_executables 2>&1 | grep -E "error" | head -30` and for each failing test file, instantiate the fixture and `attach(p)` after constructing each Particle that touches force/torque (known: `rotation_test.cpp`, `lb_particle_coupling_test.cpp`; expect a handful more, e.g. propagation or Verlet tests). For `Particle_serialization_test.cpp`: force is no longer serialized — update the expected byte layout/fields if the test asserts them.

- [ ] **Step 4: Sweep the remaining compile errors**

Run: `make -C build -j8 2>&1 | grep -E "error" | head -40`
Iterate: every error is either (a) a Task-4-pattern escape — fix with the Task-4 rules (explicit `Utils::Vector3d` conversions, no reference bindings), or (b) a genuinely new consumer (report it in your task report). `src/core/system/GpuParticleData.cpp:138-140` uses `p.force()[j] +=` — proxy `operator[]` covers it; the file only compiles with CUDA, so ALSO run a grep-review of it manually and make the same-shape edits blind (mechanical only).

- [ ] **Step 5: Full test battery + identity**

```bash
make -C build -j8 && make -C build -j8 unit_tests_executables
ctest --test-dir build -L unit_test --output-on-failure
ctest --test-dir build -R "particle|cellsystem|integrator|langevin|brownian|npt|rotation|coulomb|dipol|lees_edwards|virtual_sites|collision|lb" --output-on-failure
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode lj
./build/pypresso maintainer/benchmarks/trajectory_identity.py --mode p3m
```
Expected: all tests pass; identity hashes IDENTICAL to `.superpowers/sdd/phase1-identity-reference.txt`. A hash mismatch is a stop-the-line bug (phases 2-6 are pure storage moves — spec section 5); report BLOCKED with the diff of the first divergent quantity rather than rationalizing it.

- [ ] **Step 6: Commit**

```bash
git add -A src/core src/script_interface
git commit -m "core: force and torque live in ParticleStore columns (phase 2 flip)

The ParticleForce member leaves the Particle struct in the same commit
that makes the columns authoritative (spec section 4, single ownership).
Non-const accessors return a write-through VectorReference; const
accessors return values. Python getters fetch from the owning rank.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Phase-2 checkpoint

**Files:**
- Modify: `.superpowers/sdd/progress.md` (ledger), memory update happens outside the repo.

- [ ] **Step 1: Guard + full unit tests**

```bash
./maintainer/CI/check_cell_storage_mutations.sh
ctest --test-dir build -L unit_test --output-on-failure
```
Expected: guard OK; unit tests green.

- [ ] **Step 2: Full Python suite**

```bash
make -C build -j8 check_python_skip_long
make -C build -j8 check_python
```
Expected: both fully green. A failing test: retry once alone (statistical long tests); deterministic failure → BLOCKED with the log tail, do not weaken tests.

- [ ] **Step 3: Identity gate (final)**

Both modes, compare to `.superpowers/sdd/phase1-identity-reference.txt`: IDENTICAL.

- [ ] **Step 4: Benchmark gate** (machine must be quiet; the gate enforces foreign-CPU check; never bypass)

```bash
python3 maintainer/benchmarks/benchmark_gate.py check-load
python3 maintainer/benchmarks/benchmark_gate.py run --pypresso build/pypresso --output /tmp/phase2-current.csv --repetitions 3
python3 maintainer/benchmarks/benchmark_gate.py compare --baseline maintainer/benchmarks/baselines/phase0-baseline.csv --current /tmp/phase2-current.csv
```
Expected: PASS (≤5% cumulative vs phase-0). Watch items from phase 1: p3m 1-rank (was 1.046) and 4-rank (1.036) — report their new ratios explicitly. If the gate FAILS: report BLOCKED with the table; likely suspects are proxy overhead in `init_forces_and_thermostat` / the contribute-back loop — do not tune anything without review.

- [ ] **Step 5: Ledger**

Append to `.superpowers/sdd/progress.md`: `Phase 2: complete (commits <first>..<last>, identity OK, gate table summary)`.

---

## Plan Self-Review Notes

- **Spec coverage:** row bookkeeping + `store_index` + cell-sorted rebuild + dirty marking (Task 3); force/torque eviction with ScatterView contribution preserved (Tasks 4-7; the ScatterView machinery is untouched and writes through the flipped accessors); per-field ghost-force reduction (Task 5 — force is serialized as explicit field values, never through particle structs; full columnar gather arrives with the ghost-comm overhaul in phase 3); single-ownership commit discipline (Task 7 step 1); numerical identity gate (Tasks 1, 8; spec section 5); ≤5% benchmark budget (Task 8).
- **Deliberate scope cuts:** no unification of `m_local_force` scratch with the store column (different index spaces — `m_unique_particles` dedups ghosts; revisit in phase 3), no device-side columns beyond what DualView gives for free, no per-cell row bookkeeping (mark-dirty/full-rebuild per spec).
- **Type consistency check:** `VectorReference(double*, std::size_t)` (Tasks 2, 3); `force_reference(int)`, `force_value(int)`, `attach_to_store(ParticleStore&, int)`, `store_row()`, `ensure_particle_store_synchronized()` used identically in Tasks 3, 6, 7.
- **Line numbers** are as of commit `bb9072bf6c`; the code snippets are authoritative when lines have drifted.
