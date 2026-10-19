# ParticleStore Migration — Phases 0 & 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish the performance-gate infrastructure (phase 0) and centralize every particle insert/remove/move behind a small primitives layer (phase 1) so later phases have single hook points for mirroring rows into the `ParticleStore`.

**Architecture:** Phase 0 builds the worktree, adds a benchmark gate script with a shared-machine load check, records committed LJ/P3M baselines, and scaffolds `src/core/particle_store/`. Phase 1 is a pure refactor (zero behavior change): all direct `ParticleList` mutations on *cell storage* are rerouted through new `CellParticleStorage::` primitives; mutations of plain communication buffers stay raw. A grep-based guard script prevents regressions.

**Tech Stack:** C++20, CMake, Boost.Test (`espresso_unit_test` macro), Python 3 (benchmark gate), ctest.

**Spec:** `docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md`

## Global Constraints

- Work only in the worktree `/tikhome/weeber/es/.claude/worktrees/eliminate_Particle` (branch `worktree-eliminate_Particle`). Never touch other worktrees or `/tikhome/weeber/es` itself. All commands below run from the worktree root unless stated.
- Build with `make -j8` — never `-j$(nproc)`.
- Use `git -C <path>` forms or run git from the worktree root — do not chain `cd X && git ...`.
- Tests must be green at every commit. Success bar (from spec): full Python test suite passes; at most **5% cumulative** regression relative to the phase-0 baseline on LJ and P3M benchmarks at (1) `--particles_per_core 1000` on 1 and 4 MPI ranks, (2) `--particles_per_core 4000` with 4 OMP threads.
- **Shared machine:** never record benchmark timings while other users' heavy jobs run. The gate script enforces a load check; do not bypass it with `--force` when recording baselines or gate decisions.
- Phase 1 is a **pure refactor**: no observable behavior change, no performance-relevant change (all primitives are inline forwarding functions).
- Identifier naming: full words, no abbreviations (e.g. `insert_particle`, not `ins_part`).
- Do not re-run cmake configure once `build/` contains CMake build files (repo convention, AGENT.md), except in Task 1 where the directory is first created.
- A pre-commit hook runs clang-format on committed C++ files; if it modifies files, re-stage and commit again.
- Commit messages end with:
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`

---

### Task 1: Configure and build ESPResSo in the worktree

**Files:**
- Create: `build/` (build artifacts only — never committed; `/build*` is in `.gitignore`)

**Interfaces:**
- Produces: `build/pypresso` (Python launcher), `make -C build <target>` targets used by all later tasks: `check_python_skip_long`, `check_python`, `unit_tests_executables`.

- [ ] **Step 1: Configure (Release, defaults; CUDA stays OFF by default)**

Run: `cmake -B build -S . -D CMAKE_BUILD_TYPE=Release`
Expected: configuration succeeds, summary shows Kokkos and Cabana found. If a compiler error involving CUDA appears (it should not — `ESPRESSO_BUILD_WITH_CUDA` defaults to `OFF`), reconfigure with clang: `cmake -B build -S . -D CMAKE_BUILD_TYPE=Release -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++` (fresh `build/` directory first).

- [ ] **Step 2: Build**

Run: `make -C build -j8`
Expected: completes without errors; `build/pypresso` exists. (Takes tens of minutes.)

- [ ] **Step 3: Smoke-test the Python interface**

```bash
printf 'import espressomd\nsystem = espressomd.System(box_l=[10, 10, 10])\nsystem.part.add(pos=[1, 2, 3])\nprint("smoke-ok", len(system.part))\n' > /tmp/smoke.py
./build/pypresso /tmp/smoke.py
```
Expected output contains: `smoke-ok 1`

- [ ] **Step 4: Build and run the existing C++ unit tests**

Run: `make -C build -j8 unit_tests_executables && ctest --test-dir build -L unit_test --output-on-failure`
Expected: `100% tests passed`

- [ ] **Step 5: No commit** (build artifacts only). Verify `git -C . status --short` shows no unexpected tracked changes.

---

### Task 2: Benchmark gate script with shared-machine load check

**Files:**
- Create: `maintainer/benchmarks/benchmark_gate.py`

**Interfaces:**
- Consumes: `build/pypresso`, `maintainer/benchmarks/lj.py`, `maintainer/benchmarks/p3m.py` (both accept `--particles_per_core N --output FILE`; they append CSV rows with header `"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"` — see `maintainer/benchmarks/benchmarks.py::write_report`).
- Produces: CLI used by Tasks 3 and 10 and by every later phase checkpoint:
  - `python3 maintainer/benchmarks/benchmark_gate.py check-load [--max-load FLOAT]` → exit 0 quiet, exit 2 busy
  - `python3 maintainer/benchmarks/benchmark_gate.py run --pypresso PATH --output CSV [--repetitions N]` → runs the 6 configs × N repetitions, appends rows to CSV; exit 2 if machine busy
  - `python3 maintainer/benchmarks/benchmark_gate.py compare --baseline CSV --current CSV [--max-regression 0.05]` → exit 0 pass, exit 1 regression

- [ ] **Step 1: Write the script**

```python
#!/usr/bin/env python3
#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Benchmark gate for the ParticleStore migration.

Runs the agreed LJ/P3M benchmark configurations, records timings to CSV
(same format as maintainer/benchmarks/benchmarks.py:write_report), and
compares against a committed baseline with a regression threshold.

The machine is shared: 'check-load' refuses to measure while other jobs
are running. See docs/superpowers/specs/
2026-07-03-array-based-particle-storage-design.md, section 5.
"""

import argparse
import csv
import collections
import os
import pathlib
import subprocess
import sys

BENCHMARK_DIRECTORY = pathlib.Path(__file__).resolve().parent

# The six agreed configurations (spec: success bar).
# fields: script, extra script arguments, MPI ranks, OMP threads
CONFIGURATIONS = [
    ("lj.py", ["--particles_per_core", "1000"], 1, 1),
    ("lj.py", ["--particles_per_core", "1000"], 4, 1),
    ("lj.py", ["--particles_per_core", "4000"], 1, 4),
    ("p3m.py", ["--particles_per_core", "1000"], 1, 1),
    ("p3m.py", ["--particles_per_core", "1000"], 4, 1),
    ("p3m.py", ["--particles_per_core", "4000"], 1, 4),
]


def machine_load():
    with open("/proc/loadavg") as f:
        return float(f.read().split()[0])


def check_load(max_load):
    load = machine_load()
    if load > max_load:
        print(f"REFUSING to benchmark: 1-minute load average {load:.2f} "
              f"exceeds threshold {max_load:.2f} (shared machine).")
        print("Top CPU consumers:")
        subprocess.run(["ps", "-eo", "user,pcpu,pmem,comm", "--sort=-pcpu"],
                       check=False, stdout=sys.stdout)
        return False
    print(f"Machine quiet: load average {load:.2f} <= {max_load:.2f}.")
    return True


def run_benchmarks(pypresso, output, repetitions, max_load):
    pypresso = pathlib.Path(pypresso).resolve()
    output = pathlib.Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    for script, arguments, n_ranks, n_threads in CONFIGURATIONS:
        for repetition in range(repetitions):
            if not check_load(max_load):
                return 2
            command = ["mpirun", "-n", str(n_ranks), str(pypresso),
                       str(BENCHMARK_DIRECTORY / script)] + arguments + \
                      [f"--output={output}"]
            environment = dict(os.environ)
            environment["OMP_NUM_THREADS"] = str(n_threads)
            environment["OMP_PROC_BIND"] = "false"
            print(f"[{script} {' '.join(arguments)} ranks={n_ranks} "
                  f"threads={n_threads}] repetition {repetition + 1}"
                  f"/{repetitions}")
            subprocess.run(command, check=True, env=environment)
    print(f"Wrote timings to {output}")
    return 0


def read_rows(filepath):
    """Group rows by configuration; keep the minimum 'mean' per group."""
    groups = collections.OrderedDict()
    with open(filepath, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["script"], row["arguments"],
                   row["ranks"], row["threads"], row["label"])
            mean = float(row["mean"])
            if key not in groups or mean < groups[key]:
                groups[key] = mean
    return groups


def compare(baseline_path, current_path, max_regression):
    baseline = read_rows(baseline_path)
    current = read_rows(current_path)
    missing = [k for k in baseline if k not in current]
    if missing:
        print(f"FAIL: configurations missing from {current_path}:")
        for key in missing:
            print(f"  {key}")
        return 1
    failed = False
    print(f"{'configuration':70s} {'baseline':>12s} {'current':>12s} "
          f"{'ratio':>8s}")
    for key, baseline_mean in baseline.items():
        current_mean = current[key]
        ratio = current_mean / baseline_mean
        marker = ""
        if ratio > 1.0 + max_regression:
            marker = "  <-- REGRESSION"
            failed = True
        name = f"{key[0]} {key[1]} ranks={key[2]} threads={key[3]}"
        print(f"{name:70s} {baseline_mean:12.3e} {current_mean:12.3e} "
              f"{ratio:8.3f}{marker}")
    if failed:
        print(f"FAIL: regression beyond {100 * max_regression:.0f}% budget.")
        return 1
    print(f"PASS: all configurations within {100 * max_regression:.0f}% "
          f"of baseline.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_load = subparsers.add_parser("check-load")
    parser_load.add_argument("--max-load", type=float, default=2.0)

    parser_run = subparsers.add_parser("run")
    parser_run.add_argument("--pypresso", required=True)
    parser_run.add_argument("--output", required=True)
    parser_run.add_argument("--repetitions", type=int, default=3)
    parser_run.add_argument("--max-load", type=float, default=2.0)

    parser_compare = subparsers.add_parser("compare")
    parser_compare.add_argument("--baseline", required=True)
    parser_compare.add_argument("--current", required=True)
    parser_compare.add_argument("--max-regression", type=float, default=0.05)

    args = parser.parse_args()
    if args.command == "check-load":
        return 0 if check_load(args.max_load) else 2
    if args.command == "run":
        return run_benchmarks(args.pypresso, args.output, args.repetitions,
                              args.max_load)
    if args.command == "compare":
        return compare(args.baseline, args.current, args.max_regression)
    return 1


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Test `compare` with synthetic CSVs — pass case**

```bash
cat > /tmp/base.csv << 'EOF'
"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"
"lj.py","--particles_per_core 1000",1,1,1.000e-04,1.0e-06,5000,10.0,""
"lj.py","--particles_per_core 1000",1,1,1.100e-04,1.0e-06,5000,10.0,""
EOF
cat > /tmp/curr.csv << 'EOF'
"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"
"lj.py","--particles_per_core 1000",1,1,1.030e-04,1.0e-06,5000,10.0,""
EOF
python3 maintainer/benchmarks/benchmark_gate.py compare --baseline /tmp/base.csv --current /tmp/curr.csv; echo "exit=$?"
```
Expected: table shows ratio `1.030` (min-of-means: 1.00e-04 baseline), `PASS`, `exit=0`

- [ ] **Step 3: Test `compare` — regression case**

```bash
cat > /tmp/curr_bad.csv << 'EOF'
"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"
"lj.py","--particles_per_core 1000",1,1,1.200e-04,1.0e-06,5000,10.0,""
EOF
python3 maintainer/benchmarks/benchmark_gate.py compare --baseline /tmp/base.csv --current /tmp/curr_bad.csv; echo "exit=$?"
```
Expected: `REGRESSION` marker (ratio 1.200), `FAIL`, `exit=1`

- [ ] **Step 4: Test `check-load` refusal path**

Run: `python3 maintainer/benchmarks/benchmark_gate.py check-load --max-load 0.0001; echo "exit=$?"`
Expected: `REFUSING to benchmark`, process table, `exit=2`

- [ ] **Step 5: Test `check-load` quiet path**

Run: `python3 maintainer/benchmarks/benchmark_gate.py check-load --max-load 100; echo "exit=$?"`
Expected: `Machine quiet`, `exit=0`

- [ ] **Step 6: Commit**

```bash
git add maintainer/benchmarks/benchmark_gate.py
git commit -m "benchmarks: add benchmark gate with shared-machine load check

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Record the phase-0 benchmark baseline

**Files:**
- Create: `maintainer/benchmarks/baselines/phase0-baseline.csv`
- Create: `maintainer/benchmarks/baselines/phase0-machine-info.txt`

**Interfaces:**
- Consumes: `benchmark_gate.py run` (Task 2), `build/pypresso` (Task 1).
- Produces: the committed baseline CSV that ALL later phase checkpoints compare against (`compare --baseline maintainer/benchmarks/baselines/phase0-baseline.csv`).

- [ ] **Step 1: Verify the machine is quiet**

Run: `python3 maintainer/benchmarks/benchmark_gate.py check-load`
Expected: `Machine quiet`. If not, STOP and wait — do not record a baseline on a loaded machine.

- [ ] **Step 2: Run the baseline (6 configs × 3 repetitions; expect roughly 30–60 minutes)**

Run: `python3 maintainer/benchmarks/benchmark_gate.py run --pypresso build/pypresso --output maintainer/benchmarks/baselines/phase0-baseline.csv --repetitions 3`
Expected: exit 0; CSV contains 18 data rows (`wc -l` = 19 with header).

- [ ] **Step 3: Record machine info**

```bash
{
  echo "date: $(date -Iseconds)"
  echo "host: $(hostname)"
  echo "commit: $(git -C . rev-parse HEAD)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs)"
  echo "cores: $(nproc)"
  echo "compiler: $(grep -m1 CMAKE_CXX_COMPILER: build/CMakeCache.txt)"
  echo "build_type: $(grep -m1 CMAKE_BUILD_TYPE: build/CMakeCache.txt)"
} > maintainer/benchmarks/baselines/phase0-machine-info.txt
cat maintainer/benchmarks/baselines/phase0-machine-info.txt
```
Expected: all fields populated.

- [ ] **Step 4: Sanity-check the gate against itself**

Run: `python3 maintainer/benchmarks/benchmark_gate.py compare --baseline maintainer/benchmarks/baselines/phase0-baseline.csv --current maintainer/benchmarks/baselines/phase0-baseline.csv`
Expected: all ratios `1.000`, `PASS`.

- [ ] **Step 5: Commit**

```bash
git add maintainer/benchmarks/baselines/
git commit -m "benchmarks: record phase-0 baseline for ParticleStore migration

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: ParticleStore scaffolding

**Files:**
- Create: `src/core/particle_store/ParticleStore.hpp`
- Create: `src/core/unit_tests/ParticleStore_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt` (add one `espresso_unit_test` line next to the existing ones)

**Interfaces:**
- Produces: `class ParticleStore` (header-only for now) with `number_of_local_particles()` and `number_of_ghost_particles()`, both `std::size_t`. Phase 2 will grow this class; this task only proves the directory, include-path, and unit-test wiring.

- [ ] **Step 1: Write the failing test**

`src/core/unit_tests/ParticleStore_test.cpp` (copy the GPL header block from `src/core/unit_tests/Particle_test.cpp` lines 1–18, then):

```cpp
#define BOOST_TEST_MODULE ParticleStore test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "particle_store/ParticleStore.hpp"

BOOST_AUTO_TEST_CASE(default_constructed_store_is_empty) {
  ParticleStore const store{};
  BOOST_CHECK_EQUAL(store.number_of_local_particles(), 0ul);
  BOOST_CHECK_EQUAL(store.number_of_ghost_particles(), 0ul);
}
```

Register in `src/core/unit_tests/CMakeLists.txt` (alphabetically near the `Particle_test.cpp` entry, matching its `DEPENDS` clause — check the existing line and copy its dependencies):

```cmake
espresso_unit_test(SRC ParticleStore_test.cpp DEPENDS espresso::core)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `make -C build -j8 ParticleStore_test 2>&1 | tail -5`
Expected: FAIL — `particle_store/ParticleStore.hpp: No such file or directory`
(If the target does not exist yet, cmake re-runs automatically on `make`; if not, run `cmake -B build -S .` once to pick up the new test file.)

- [ ] **Step 3: Write the header**

`src/core/particle_store/ParticleStore.hpp` (copy the GPL header block from `src/core/Particle.hpp` lines 1–18, then):

```cpp
#pragma once

#include <cstddef>

/**
 * @brief Array-based particle storage (structure of arrays).
 *
 * Owns every per-particle quantity in a single index space: local
 * particles first, contiguously sorted by cell, ghost particles appended
 * after the locals. Fields are grouped into parameters, state, and
 * observables. See
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 *
 * This is scaffolding (migration phase 0): the columns are added in
 * migration phase 2.
 */
class ParticleStore {
public:
  auto number_of_local_particles() const { return m_number_of_local_particles; }
  auto number_of_ghost_particles() const { return m_number_of_ghost_particles; }

private:
  std::size_t m_number_of_local_particles = 0u;
  std::size_t m_number_of_ghost_particles = 0u;
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `make -C build -j8 ParticleStore_test && ctest --test-dir build -R ParticleStore_test --output-on-failure`
Expected: PASS (`100% tests passed`)

- [ ] **Step 5: Commit**

```bash
git add src/core/particle_store/ParticleStore.hpp src/core/unit_tests/ParticleStore_test.cpp src/core/unit_tests/CMakeLists.txt
git commit -m "core: scaffold ParticleStore (migration phase 0)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: CellParticleStorage primitives (TDD)

**Files:**
- Create: `src/core/cell_system/ParticleListOperations.hpp`
- Create: `src/core/unit_tests/ParticleListOperations_test.cpp`
- Modify: `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `ParticleList` = `Utils::Bag<Particle>` (`src/core/ParticleList.hpp`); `Utils::Bag` API: `T &insert(T &&)`, `iterator erase(iterator)` (swap-with-back, order NOT preserved), `void clear()`, `void resize(std::size_t)` (`src/utils/include/utils/Bag.hpp:118-171`).
- Produces (used verbatim by Tasks 6–9): namespace `CellParticleStorage` with
  - `Particle &insert_particle(ParticleList &storage, Particle &&particle)`
  - `std::pair<Particle, ParticleList::iterator> extract_particle(ParticleList &storage, ParticleList::iterator position)`
  - `ParticleList::iterator erase_particle(ParticleList &storage, ParticleList::iterator position)`
  - `void clear_particles(ParticleList &storage)`
  - `void resize_ghost_storage(ParticleList &storage, std::size_t count)`

- [ ] **Step 1: Write the failing test**

`src/core/unit_tests/ParticleListOperations_test.cpp` (GPL header block as in Task 4, then):

```cpp
#define BOOST_TEST_MODULE ParticleListOperations test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "ParticleList.hpp"
#include "cell_system/ParticleListOperations.hpp"

#include <utility>

static Particle make_particle(int id) {
  Particle p{};
  p.id() = id;
  return p;
}

BOOST_AUTO_TEST_CASE(insert_particle_appends_and_returns_reference) {
  ParticleList storage;
  auto &stored = CellParticleStorage::insert_particle(storage, make_particle(7));
  BOOST_CHECK_EQUAL(storage.size(), 1ul);
  BOOST_CHECK_EQUAL(stored.id(), 7);
}

BOOST_AUTO_TEST_CASE(extract_particle_moves_out_with_swap_with_back) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::insert_particle(storage, make_particle(2));
  CellParticleStorage::insert_particle(storage, make_particle(3));

  auto [extracted, next] =
      CellParticleStorage::extract_particle(storage, storage.begin());
  BOOST_CHECK_EQUAL(extracted.id(), 1);
  BOOST_CHECK_EQUAL(storage.size(), 2ul);
  // swap-with-back: the last element (id 3) now occupies the freed slot
  BOOST_CHECK_EQUAL(next->id(), 3);
}

BOOST_AUTO_TEST_CASE(erase_particle_destroys_element) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::insert_particle(storage, make_particle(2));
  auto it = CellParticleStorage::erase_particle(storage, storage.begin());
  BOOST_CHECK_EQUAL(storage.size(), 1ul);
  BOOST_CHECK_EQUAL(it->id(), 2);
}

BOOST_AUTO_TEST_CASE(clear_particles_empties_storage) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::clear_particles(storage);
  BOOST_CHECK_EQUAL(storage.size(), 0ul);
}

BOOST_AUTO_TEST_CASE(resize_ghost_storage_marks_all_particles_as_ghosts) {
  ParticleList storage;
  CellParticleStorage::insert_particle(storage, make_particle(1));
  CellParticleStorage::resize_ghost_storage(storage, 3ul);
  BOOST_CHECK_EQUAL(storage.size(), 3ul);
  for (auto const &p : storage) {
    BOOST_CHECK(p.is_ghost());
  }
}
```

Register in `src/core/unit_tests/CMakeLists.txt` (same `DEPENDS` as the `Particle_test.cpp` entry):

```cmake
espresso_unit_test(SRC ParticleListOperations_test.cpp DEPENDS espresso::core)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `make -C build -j8 ParticleListOperations_test 2>&1 | tail -5`
Expected: FAIL — `cell_system/ParticleListOperations.hpp: No such file or directory`

- [ ] **Step 3: Write the header**

`src/core/cell_system/ParticleListOperations.hpp` (GPL header block, then):

```cpp
#pragma once

#include "Particle.hpp"
#include "ParticleList.hpp"

#include <cstddef>
#include <utility>

/**
 * @brief Primitives for mutating the particle storage of cells.
 *
 * Every insertion, removal, or move of a @ref Particle in *cell storage*
 * (a @ref ParticleList owned by a @ref Cell, including ghost cells) MUST
 * go through these functions. They are the single hook point for
 * mirroring rows into the ParticleStore in later migration phases; see
 * docs/superpowers/specs/2026-07-03-array-based-particle-storage-design.md
 * section 3, phase 1.
 *
 * Plain communication buffers (send/receive @ref ParticleList instances
 * that do not belong to a cell) are exempt and use the @ref Utils::Bag
 * API directly.
 *
 * maintainer/CI/check_cell_storage_mutations.sh enforces this rule.
 */
namespace CellParticleStorage {

/**
 * @brief Insert a particle into a cell's particle storage.
 * May reallocate the storage, invalidating pointers and iterators into it.
 * @return Reference to the stored particle.
 */
inline Particle &insert_particle(ParticleList &storage, Particle &&particle) {
  return storage.insert(std::move(particle));
}

/**
 * @brief Move a particle out of a cell's particle storage.
 * Swap-with-back removal: element order is not preserved.
 * @return The extracted particle and the iterator past the removed element.
 */
inline std::pair<Particle, ParticleList::iterator>
extract_particle(ParticleList &storage, ParticleList::iterator position) {
  auto particle = std::move(*position);
  auto next = storage.erase(position);
  return {std::move(particle), next};
}

/**
 * @brief Erase a particle from a cell's particle storage, destroying it.
 * Swap-with-back removal: element order is not preserved.
 * @return Iterator past the removed element.
 */
inline ParticleList::iterator erase_particle(ParticleList &storage,
                                             ParticleList::iterator position) {
  return storage.erase(position);
}

/** @brief Remove all particles from a cell's particle storage. */
inline void clear_particles(ParticleList &storage) { storage.clear(); }

/**
 * @brief Resize ghost-cell storage to exactly @p count particles.
 * Newly created particles are default-constructed; all particles in the
 * storage are (re)marked as ghosts.
 */
inline void resize_ghost_storage(ParticleList &storage, std::size_t count) {
  storage.resize(count);
  for (auto &particle : storage) {
    particle.set_ghost(true);
  }
}

} // namespace CellParticleStorage
```

Note: if `ParticleList::iterator` does not compile, check `Utils::Bag`'s iterator typedef name in `src/utils/include/utils/Bag.hpp` and use `Utils::Bag<Particle>::iterator` spelled the way `Bag` exposes it.

- [ ] **Step 4: Run test to verify it passes**

Run: `make -C build -j8 ParticleListOperations_test && ctest --test-dir build -R ParticleListOperations_test --output-on-failure`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/core/cell_system/ParticleListOperations.hpp src/core/unit_tests/ParticleListOperations_test.cpp src/core/unit_tests/CMakeLists.txt
git commit -m "core: add CellParticleStorage primitives (migration phase 1)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Reroute CellStructure through the primitives

**Files:**
- Modify: `src/core/cell_system/CellStructure.hpp:273-287` (`append_indexed_particle`)
- Modify: `src/core/cell_system/CellStructure.cpp:372-396` (`remove_particle`), `:451-458` (`remove_all_particles`)

**Interfaces:**
- Consumes: `CellParticleStorage::insert_particle`, `erase_particle`, `clear_particles` (Task 5).
- Produces: no API change — pure internal reroute.

- [ ] **Step 1: Edit `CellStructure.hpp`**

Add `#include "cell_system/ParticleListOperations.hpp"` to the include block of `CellStructure.hpp`. In `append_indexed_particle` replace

```cpp
    auto &new_part = pl.insert(std::move(p));
```

with

```cpp
    auto &new_part = CellParticleStorage::insert_particle(pl, std::move(p));
```

- [ ] **Step 2: Edit `CellStructure.cpp`**

In `remove_particle` (line ~387) replace

```cpp
        it = parts.erase(it);
```

with

```cpp
        it = CellParticleStorage::erase_particle(parts, it);
```

In `remove_all_particles` (line ~453) replace

```cpp
    cell->particles().clear();
```

with

```cpp
    CellParticleStorage::clear_particles(cell->particles());
```

(`CellStructure.cpp` includes `CellStructure.hpp`, which now provides the header.)

- [ ] **Step 3: Rebuild and run C++ unit tests**

Run: `make -C build -j8 && make -C build -j8 unit_tests_executables && ctest --test-dir build -L unit_test --output-on-failure`
Expected: build clean, `100% tests passed`

- [ ] **Step 4: Run targeted Python tests**

Run: `make -C build -j8 pypresso python_test_data && ctest --test-dir build -R "particle" --output-on-failure`
Expected: all matched tests pass (covers add/remove particle paths, serial and MPI variants).

- [ ] **Step 5: Commit**

```bash
git add src/core/cell_system/CellStructure.hpp src/core/cell_system/CellStructure.cpp
git commit -m "core: route CellStructure storage mutations through CellParticleStorage

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Reroute RegularDecomposition

**Files:**
- Modify: `src/core/cell_system/RegularDecomposition.cpp:84-99` (`move_if_local`), `:176-240` (`resort`)

**Interfaces:**
- Consumes: `CellParticleStorage::insert_particle`, `extract_particle` (Task 5).
- Produces: no API change. `move_left_or_right` (lines 101–123) and all `send_buf_*`/`recv_buf_*`/`displaced_parts` operations intentionally stay raw `Bag` calls — they mutate communication buffers, not cell storage.

- [ ] **Step 1: Edit `move_if_local`**

Add `#include "cell_system/ParticleListOperations.hpp"` to the include block of `RegularDecomposition.cpp`. Replace (line ~91)

```cpp
      target_cell->particles().insert(std::move(part));
```

with

```cpp
      CellParticleStorage::insert_particle(target_cell->particles(),
                                           std::move(part));
```

(`rest.insert(...)` and `src.clear()` stay raw: `src` and `rest` are buffers.)

- [ ] **Step 2: Edit `resort`**

Replace (lines ~192-194)

```cpp
      auto p = std::move(*it);
      it = c->particles().erase(it);
      diff.emplace_back(ModifiedList{c->particles()});
```

with

```cpp
      auto [p, next] = CellParticleStorage::extract_particle(c->particles(), it);
      it = next;
      diff.emplace_back(ModifiedList{c->particles()});
```

Replace (line ~203)

```cpp
        target_cell->particles().insert(std::move(p));
```

with

```cpp
        CellParticleStorage::insert_particle(target_cell->particles(),
                                             std::move(p));
```

Replace (line ~235)

```cpp
      sort_cell->particles().insert(std::move(part));
```

with

```cpp
      CellParticleStorage::insert_particle(sort_cell->particles(),
                                           std::move(part));
```

(`displaced_parts.insert(std::move(p));` at line ~199 stays raw: buffer.)

- [ ] **Step 3: Rebuild and run C++ unit tests**

Run: `make -C build -j8 && ctest --test-dir build -L unit_test --output-on-failure`
Expected: build clean, all pass.

- [ ] **Step 4: Run targeted Python tests (sorting/exchange heavy)**

Run: `ctest --test-dir build -R "cellsystem|particle|integrator" --output-on-failure`
Expected: all matched tests pass (includes 4-rank MPI variants that exercise `exchange_neighbors` and `resort`).

- [ ] **Step 5: Commit**

```bash
git add src/core/cell_system/RegularDecomposition.cpp
git commit -m "core: route RegularDecomposition cell mutations through CellParticleStorage

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Reroute AtomDecomposition and HybridDecomposition

**Files:**
- Modify: `src/core/cell_system/AtomDecomposition.cpp:106-145` (`resort`)
- Modify: `src/core/cell_system/HybridDecomposition.cpp:95-161` (`resort`)

**Interfaces:**
- Consumes: `CellParticleStorage::insert_particle`, `extract_particle` (Task 5).
- Produces: no API change. `send_buf`/`recv_buf` (`std::vector<Particle>`) stay raw — buffers.

- [ ] **Step 1: Edit `AtomDecomposition.cpp`**

Add `#include "cell_system/ParticleListOperations.hpp"`. Replace (lines ~125-127)

```cpp
      diff.emplace_back(RemovedParticle{it->id()});
      send_buf.at(target_node).emplace_back(std::move(*it));
      it = local().particles().erase(it);
```

with

```cpp
      diff.emplace_back(RemovedParticle{it->id()});
      auto [p, next] =
          CellParticleStorage::extract_particle(local().particles(), it);
      send_buf.at(target_node).emplace_back(std::move(p));
      it = next;
```

Replace (line ~142)

```cpp
      local().particles().insert(std::move(p));
```

with

```cpp
      CellParticleStorage::insert_particle(local().particles(), std::move(p));
```

Watch the name collision: the received-particle loop already uses `p` (`for (auto &p : parts)`); rename the structured-binding variable to `extracted` if the compiler flags shadowing:

```cpp
      auto [extracted, next] =
          CellParticleStorage::extract_particle(local().particles(), it);
      send_buf.at(target_node).emplace_back(std::move(extracted));
      it = next;
```

- [ ] **Step 2: Edit `HybridDecomposition.cpp`**

Add `#include "cell_system/ParticleListOperations.hpp"`. Apply the same transformation four times:

Lines ~110-111 (regular → n_square):

```cpp
      auto [p, next] =
          CellParticleStorage::extract_particle(cell_rd->particles(), it);
      it = next;
```

Line ~117:

```cpp
      CellParticleStorage::insert_particle(first_local_cell->particles(),
                                           std::move(p));
```

Lines ~132-133 (n_square → regular):

```cpp
        auto [p, next] =
            CellParticleStorage::extract_particle(cell_ns->particles(), it);
        it = next;
```

Lines ~141 and ~147:

```cpp
          CellParticleStorage::insert_particle(target_cell->particles(),
                                               std::move(p));
```

```cpp
          CellParticleStorage::insert_particle(first_local_cell->particles(),
                                               std::move(p));
```

- [ ] **Step 3: Rebuild and run C++ unit tests**

Run: `make -C build -j8 && ctest --test-dir build -L unit_test --output-on-failure`
Expected: build clean, all pass.

- [ ] **Step 4: Run targeted Python tests**

Run: `ctest --test-dir build -R "hybrid|nsquare|n_square|cellsystem|particle" --output-on-failure`
Expected: all matched tests pass. (If no test matches `hybrid|nsquare|n_square`, list names with `ctest --test-dir build -N | grep -i -E "hybrid|square"` and run what exists.)

- [ ] **Step 5: Commit**

```bash
git add src/core/cell_system/AtomDecomposition.cpp src/core/cell_system/HybridDecomposition.cpp
git commit -m "core: route Atom/HybridDecomposition cell mutations through CellParticleStorage

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Reroute ghost-cell resizing

**Files:**
- Modify: `src/core/ghosts.cpp:290-322` (`prepare_ghost_cell`, `put_recv_buffer`)

**Interfaces:**
- Consumes: `CellParticleStorage::resize_ghost_storage` (Task 5).
- Produces: no API change; the static helper `prepare_ghost_cell` is deleted.

- [ ] **Step 1: Edit `ghosts.cpp`**

Add `#include "cell_system/ParticleListOperations.hpp"` to the include block. Delete the whole `prepare_ghost_cell` function (lines ~290-298):

```cpp
static void prepare_ghost_cell(ParticleList *cell, std::size_t size) {
  /* Adapt size */
  cell->resize(size);

  /* Mark particles as ghosts */
  for (auto &p : *cell) {
    p.set_ghost(true);
  }
}
```

In `put_recv_buffer` replace (line ~321)

```cpp
      prepare_ghost_cell(part_list, np);
```

with

```cpp
      CellParticleStorage::resize_ghost_storage(*part_list, np);
```

- [ ] **Step 2: Rebuild and run C++ unit tests**

Run: `make -C build -j8 && ctest --test-dir build -L unit_test --output-on-failure`
Expected: build clean, all pass.

- [ ] **Step 3: Run targeted Python tests (ghost exchange heavy)**

Run: `ctest --test-dir build -R "particle|cellsystem|lees_edwards|coulomb" --output-on-failure`
Expected: all matched tests pass (MPI variants exercise ghost communication).

- [ ] **Step 4: Commit**

```bash
git add src/core/ghosts.cpp
git commit -m "core: route ghost-cell resizing through CellParticleStorage

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Mutation guard script + phase-1 checkpoint

**Files:**
- Create: `maintainer/CI/check_cell_storage_mutations.sh`

**Interfaces:**
- Consumes: rerouted sources (Tasks 6–9), `benchmark_gate.py` + baseline (Tasks 2–3).
- Produces: guard script run at every later phase checkpoint; the completed phase-1 checkpoint (tests + benchmarks).

- [ ] **Step 1: Write the guard script**

`maintainer/CI/check_cell_storage_mutations.sh`:

```bash
#!/usr/bin/env bash
# Guard for the ParticleStore migration (see docs/superpowers/specs/
# 2026-07-03-array-based-particle-storage-design.md, section 3, phase 1):
# cell particle storage may only be mutated through the primitives in
# src/core/cell_system/ParticleListOperations.hpp. This script fails if a
# direct mutation pattern reappears in src/core (unit tests excluded:
# they construct fixtures directly).
set -u
cd "$(dirname "$0")/../.."
matches=$(grep -rn --include='*.cpp' --include='*.hpp' \
    -e 'particles()\.insert(' \
    -e 'particles()\.erase(' \
    -e 'particles()\.clear(' \
    -e 'particles()\.resize(' \
    src/core \
    | grep -v 'src/core/unit_tests/' \
    | grep -v 'src/core/cell_system/ParticleListOperations.hpp')
if [ -n "${matches}" ]; then
    echo "ERROR: direct cell-storage mutation outside CellParticleStorage:"
    echo "${matches}"
    exit 1
fi
echo "OK: no direct cell-storage mutations found."
exit 0
```

Run: `chmod +x maintainer/CI/check_cell_storage_mutations.sh`

- [ ] **Step 2: Run the guard — must pass**

Run: `./maintainer/CI/check_cell_storage_mutations.sh; echo "exit=$?"`
Expected: `OK: no direct cell-storage mutations found.`, `exit=0`.
If it reports matches, those are missed reroute sites: fix them the same way as in Tasks 6–9 (cell storage → primitives; genuine buffers → adjust the grep exclusion ONLY if the match is provably a buffer, and document why in the script).

- [ ] **Step 3: Verify the guard actually catches violations**

```bash
sed -i 's|CellParticleStorage::clear_particles(cell->particles());|cell->particles().clear();|' src/core/cell_system/CellStructure.cpp
./maintainer/CI/check_cell_storage_mutations.sh; echo "exit=$?"
git checkout -- src/core/cell_system/CellStructure.cpp
```
Expected: `ERROR: direct cell-storage mutation ...`, `exit=1`; then the revert restores the reroute.

- [ ] **Step 4: Phase-1 checkpoint — full Python test suite**

Run: `make -C build -j8 check_python_skip_long`
Expected: all tests pass (takes a while). Then the complete suite including long-labelled tests:
Run: `make -C build -j8 check_python`
Expected: all tests pass. If time pressure requires deferring the `long`-labelled statistical tests, say so explicitly in the checkpoint report — do not silently skip.

- [ ] **Step 5: Phase-1 checkpoint — benchmark gate**

```bash
python3 maintainer/benchmarks/benchmark_gate.py check-load
python3 maintainer/benchmarks/benchmark_gate.py run --pypresso build/pypresso --output /tmp/phase1-current.csv --repetitions 3
python3 maintainer/benchmarks/benchmark_gate.py compare --baseline maintainer/benchmarks/baselines/phase0-baseline.csv --current /tmp/phase1-current.csv
```
Expected: `PASS` — phase 1 is inline forwarding only, ratios should be ~1.00.

- [ ] **Step 6: Commit**

```bash
git add maintainer/CI/check_cell_storage_mutations.sh
git commit -m "CI: guard against direct cell-storage mutations (migration phase 1)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Plan Self-Review Notes

- **Spec coverage:** phase 0 = baselines (Tasks 1–3), comparison script with load check (Task 2), `particle_store` scaffolding (Task 4). Phase 1 = centralized primitives (Task 5), all seven inventoried direct-mutation groups rerouted (Tasks 6–9: CellStructure ×3, RegularDecomposition ×2 groups, AtomDecomposition, HybridDecomposition, ghosts), enforcement + checkpoint (Task 10). Buffer mutations (`move_left_or_right`, send/recv buffers, `displaced_parts`) are deliberately exempt per the primitives' contract.
- **Line numbers** are as of commit `c920cee7f3`; verify context with the shown code snippets before editing — the snippet is authoritative, not the line number.
- **Naming:** `CellParticleStorage` namespace, full-word function names, consistent across Tasks 5–10.
