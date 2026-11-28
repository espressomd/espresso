# Distributed Voxel-Mask Shape Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `VoxelMask` shape whose solid/pore mask is read per rank from HDF5 and stored distributed, usable for shape-based constraints and for LB boundaries without any global raster.

**Architecture:** A pure `Shapes::VoxelMask` (no MPI, no HDF5) stores a 1-bit voxel block for the local box plus halo and answers distance queries. A core HDF5 reader (`IO::VoxelMaskFile`, HighFive, MPI-IO driver) reads hyperslabs. The script-interface `VoxelMask` glues both, using a loader callback that the constraint system calls on node-grid and periodicity changes. A new rank-local LB boundary path evaluates `is_inside()` per local LB node for any shape.

**Tech Stack:** C++20, Boost.MPI, HighFive / parallel HDF5, waLBerla bridge, Python + h5py for tests.

**Spec:** `docs/superpowers/specs/2026-10-24-distributed-voxel-mask-shape-design.md`

## Global Constraints

- No code path materializes the full mask, or a full-grid LB raster, on one rank or in Python.
- No voxel data passes through the script interface or `GlobalContext`.
- `VoxelMask` exists only with `ESPRESSO_BUILD_WITH_HDF5` (C++ macro `ESPRESSO_H5MD`, Python feature `H5MD`).
- HDF5 dataset: 3D, indexed `[x][y][z]` in C order; nonzero means solid; mask extent = dataset shape × `voxel_size`.
- In memory: 1 bit per voxel. No extra per-voxel data.
- Periodic halo wraps by the box length. Outside the mask, and beyond a non-periodic box face, is pore. Mask voxels beyond the box are ignored. None of these is an error.
- `calculate_dist()`: signed distance to the nearest solid voxel cube within the halo; positive in pore, negative in solid; `vec = pos - closest point on the interface`.
- New code runs with 1, 2, 3 and 4 MPI ranks.
- Run `maintainer/format` on changed files before each commit.
- Build with `make -j8`, never `-j$(nproc)`.
- Out of scope: load balancing, `PotentialField`/`ForceField`, EK boundaries, `get_boundary_force_from_shape()` in the local path, non-uniform slip velocity in the local path.

## Plan-level decisions (not in the spec)

- **Independent, not collective, reads.** Each rank opens the file with the MPI-IO driver and reads its own hyperslabs with independent transfer. With a periodic wrap a rank needs up to 3 segments per axis, so ranks issue different numbers of reads, which collective transfer does not allow.
- **No automatic halo-vs-cutoff check.** The halo is documented as: `halo * voxel_size >= max_cutoff + skin / 2`. A query whose own voxel lies outside the stored region raises an error. The spec marked this check as a suggestion.
- **Local LB path replaces the global path for uniform velocity, for all shapes.** The global path stays for array-valued velocity. `VoxelMask` with array velocity raises `ValueError`.
- **`Union` containing `VoxelMask`:** no special handling. Out-of-region errors surface.

## Review Focus

1. **Particle drifted across a periodic face before resort.** Its folded position lies on the far side of the box, and the shape must use the periodic image inside the stored region. Test: Task 1, `image_position_across_periodic_face`.
2. **Missing file or dataset.** All ranks must raise the same Python exception without deadlock. Test: Task 3, `test_missing_file_raises`.
3. **Dataset written by h5py as `bool` or `int32`.** It must load the same as `uint8`. Test: Task 2, `reads_bool_and_int32_datasets`.
4. **Mask smaller or larger than the box in a periodic direction.** The uncovered part is pore, and extra voxels are ignored. Test: Task 3, `test_mask_smaller_than_box` and `test_mask_larger_than_box`.
5. **No solid within the halo.** The distance exceeds `halo * voxel_size`, and the constraint force is exactly zero. Test: Task 1, `far_from_solid`; Task 4, `test_force_zero_far_from_wall`.

---

### Task 1: Build setup and core `Shapes::VoxelMask`

**Files:**
- Create: `src/shapes/include/shapes/VoxelMask.hpp`, `src/shapes/src/VoxelMask.cpp`
- Modify: `src/shapes/include/shapes/Shape.hpp` (add virtual hook), `src/shapes/CMakeLists.txt` (add source)
- Test: `src/shapes/unit_tests/VoxelMask_test.cpp`, `src/shapes/unit_tests/CMakeLists.txt`

**Interfaces:**
- Produces, in `Shapes::Shape`: `virtual void on_local_geometry_change() {}`.
- Produces, in namespace `Shapes`:

```cpp
class VoxelMask : public Shape {
public:
  struct Geometry {
    Utils::Vector3i mask_shape;   // dataset shape, voxels
    double voxel_size;
    Utils::Vector3d origin;       // position of the mask's lower corner
    Utils::Vector3d box_l;
    std::array<bool, 3> periodic;
    int halo;                     // voxels, >= 1
  };
  struct Region { Utils::Vector3i lower; Utils::Vector3i extent; };

  // Unfolded voxel index j <-> voxel center origin + (j + 0.5) * voxel_size.
  static Region stored_region(Geometry const &g, Utils::Vector3d const &local_lower,
                              Utils::Vector3d const &local_upper);
  static int mask_index(Geometry const &g, int axis, int j); // -1 = pore

  void set_geometry(Geometry const &g);
  Geometry const &geometry() const;
  // voxels: row-major over extent, x slowest; nonzero = solid
  void set_local_data(Region const &region, std::vector<std::uint8_t> const &voxels);
  void set_loader(std::function<void(VoxelMask &)> loader);
  void on_local_geometry_change() override;   // calls loader if set

  void calculate_dist(Utils::Vector3d const &pos, double &dist,
                      Utils::Vector3d &vec) const override;
  bool is_inside(Utils::Vector3d const &pos) const override;
  std::size_t local_voxel_count() const;
  std::size_t local_memory_bytes() const;
};
```

Exact rules (the signatures do not determine these):

```text
stored_region: lower[a] = floor((local_lower[a] - origin[a]) / h) - halo
               upper[a] = ceil ((local_upper[a] - origin[a]) / h) + halo   (exclusive)
mask_index:    c = origin[a] + (j + 0.5) h
               periodic[a]: c -= floor(c / box_l[a]) * box_l[a]
               else if c < 0 or c >= box_l[a]: return -1
               m = floor((c - origin[a]) / h); return 0 <= m < mask_shape[a] ? m : -1
image choice:  for each periodic axis try pos[a] + n*box_l[a], n in {0, -1, +1};
               take the first whose voxel lies in the stored region;
               none -> std::runtime_error("VoxelMask: position outside the local region of this rank")
search:        voxels within +-halo of the query voxel, clipped to the stored region
pore query:    min over solid cubes of |pos - clamp(pos, cube)|
solid query:   -(min over pore cubes of |pos - clamp(pos, cube)|)
none found:    dist = (halo + 1) * h, vec = dist * (1, 0, 0)
is_inside:     the voxel containing pos is solid
```

- [ ] **Step 1: Configure a build in the worktree**

Run: `cmake -S . -B build -DESPRESSO_BUILD_WITH_HDF5=ON -DESPRESSO_BUILD_WITH_WALBERLA=ON -DESPRESSO_BUILD_WITH_CUDA=OFF -DESPRESSO_BUILD_TESTS=ON && make -C build -j8`
Expected: build succeeds. This build is reused by all tasks.

- [ ] **Step 2: Write the failing unit tests** in `VoxelMask_test.cpp` (Boost.Test, style of `Wall_test.cpp`). Use h = 1, origin = 0, box_l = 10, mask_shape = 10³, halo = 2, all axes periodic, and a local box of [0, 5)³ unless stated.
  - `stored_region_bounds`: lower = (-2,-2,-2), extent = (9,9,9).
  - `mask_index_periodic_wrap`: `mask_index(g, 0, -1) == 9`, `mask_index(g, 0, 10) == 0`; with `periodic[0] = false`: `-1` for both.
  - `mask_index_box_not_multiple_of_voxel`: box_l = 10.5, mask 11 voxels: `mask_index(g, 0, 9) == 9`, `mask_index(g, 0, 10) == 0` (center 10.5 folds to 0), `mask_index(g, 0, 11) == 1`.
  - `slab_distance`: solid where mask x < 2. At pos (3.5, 2, 2): dist = 1.5, vec = (1.5, 0, 0). At pos (1.5, 2, 2): dist = -0.5, vec = (-0.5, 0, 0). `is_inside` matches the sign.
  - `edge_distance`: a single solid voxel at (2,2,2); pos (3.5, 3.5, 2.5): dist = √0.5, vec = (0.5, 0.5, 0).
  - `far_from_solid`: an empty mask: dist = 3.0.
  - `image_position_across_periodic_face`: the plane mask x = 9 is solid, everything else pore; local box [0,5)³. Folded pos (9.8, 2.5, 2.5) lies outside the stored region; its image (-0.2, 2.5, 2.5) lies in stored voxel -1, which maps to mask x = 9. Expect `is_inside` true, dist = -0.2, vec = (-0.2, 0, 0). Pos (0.3, 2.5, 2.5): dist = 0.3, vec = (0.3, 0, 0).
  - `outside_region_throws`: pos (7.5, 2, 2) with no periodic image in range: `BOOST_CHECK_THROW(..., std::runtime_error)`.
  - `memory_is_bit_packed`: `local_memory_bytes() <= local_voxel_count() / 8 + 64`.
  - `loader_called_on_geometry_change`: a loader that counts calls; `on_local_geometry_change()` increments it.

- [ ] **Step 3: Register and run; verify failure.** Add `espresso_unit_test(SRC VoxelMask_test.cpp DEPENDS espresso::shapes espresso::utils)`.
Run: `make -C build -j8 VoxelMask_test`
Expected: compile failure, because `shapes/VoxelMask.hpp` does not exist.

- [ ] **Step 4: Implement `Shapes::VoxelMask`** using the rules above. Store bits in `std::vector<std::uint64_t>`.

- [ ] **Step 5: Run the tests.**
Run: `make -C build -j8 VoxelMask_test && ctest --test-dir build -R VoxelMask_test --output-on-failure`
Expected: PASS.

- [ ] **Step 6: Commit.** `git add src/shapes && git commit -m "shapes: Add VoxelMask shape core"`

---

### Task 2: HDF5 mask reader

**Files:**
- Create: `src/core/io/voxel_mask/VoxelMaskFile.hpp`, `src/core/io/voxel_mask/VoxelMaskFile.cpp`, `src/core/io/voxel_mask/CMakeLists.txt`
- Modify: `src/core/io/CMakeLists.txt` (`add_subdirectory(voxel_mask)` after `writer`; append the source to `espresso_hdf5` with `target_sources` inside `if(ESPRESSO_BUILD_WITH_HDF5)`)
- Test: `src/core/unit_tests/VoxelMaskFile_test.cpp`, `src/core/unit_tests/CMakeLists.txt`

**Interfaces:**
- Produces, in namespace `IO`:

```cpp
class VoxelMaskFile {
public:
  // Opens read-only with the MPI-IO driver; throws std::runtime_error naming
  // the path or dataset when either is missing or the dataset is not 3D.
  VoxelMaskFile(std::string const &path, std::string const &dataset,
                boost::mpi::communicator const &comm);
  Utils::Vector3i shape() const;
  // Independent read of [offset, offset + count); row-major, x slowest;
  // values converted to 0/1 (nonzero -> 1).
  std::vector<std::uint8_t> read(Utils::Vector3i const &offset,
                                 Utils::Vector3i const &count) const;
};
```

- [ ] **Step 1: Write the failing test** (Boost.Test, `NUM_PROC 2`). Rank 0 writes `$TMPDIR`-relative test files with HighFive, followed by a barrier.
  - `reads_hyperslab`: a 4×5×6 `uint8` dataset with value `(x*100 + y*10 + z) % 2`. Each rank reads offset (1,2,3), count (2,2,2) and checks every element against the formula.
  - `reads_bool_and_int32_datasets`: datasets of type `int32` (values 0/7) and HDF5 enum bool (as h5py writes it, via `HighFive::create_enum_boolean`). Both read as 0/1.
  - `missing_dataset_throws`: `BOOST_CHECK_THROW(IO::VoxelMaskFile(path, "nope", comm), std::runtime_error)`.

- [ ] **Step 2: Register and verify failure.** `espresso_unit_test(SRC VoxelMaskFile_test.cpp DEPENDS espresso::core Boost::mpi MPI::MPI_CXX NUM_PROC 2)` inside `if(ESPRESSO_BUILD_WITH_HDF5)`.
Run: `make -C build -j8 VoxelMaskFile_test`
Expected: compile failure, header missing.

- [ ] **Step 3: Implement `IO::VoxelMaskFile`** with HighFive (`FileAccessProps` + `MPIOFileAccess`, `DataSet::select(offset, count)`). Read into a native integer type and convert. For enum bool, read as `int8`.

- [ ] **Step 4: Run.**
Run: `make -C build -j8 VoxelMaskFile_test && ctest --test-dir build -R VoxelMaskFile_test --output-on-failure`
Expected: PASS.

- [ ] **Step 5: Commit.** `git commit -m "io: Add parallel HDF5 voxel-mask reader"`

---

### Task 3: Script interface, Python class, distributed queries

**Files:**
- Create: `src/script_interface/shapes/VoxelMask.hpp`
- Modify: `src/script_interface/shapes/initialize.cpp` (register `"Shapes::VoxelMask"` under `#ifdef ESPRESSO_H5MD`), `src/python/espressomd/shapes.py`, `doc/sphinx/constraints.rst` (a `VoxelMask` subsection after `Wall`, including the halo rule and the HDF5 layout)
- Test: `testsuite/python/shape_voxel_mask.py`, `testsuite/python/CMakeLists.txt`

**Interfaces:**
- Consumes: `Shapes::VoxelMask` (Task 1), `IO::VoxelMaskFile` (Task 2).
- Produces: script object `Shapes::VoxelMask` with read-only parameters `file` (str), `dataset` (str), `voxel_size` (float), `origin` (3 floats), `halo` (int), `mask_shape` (3 ints, from the file). Methods:
  - `calc_distance(position)`, `is_inside(position)`: the rank whose local box contains the folded position evaluates, and the result is returned on the head rank.
  - `rasterize`: raises `RuntimeError("VoxelMask cannot be rasterized globally")`.
  - `get_local_voxel_counts()`, `get_local_memory_bytes()`: one entry per rank, in rank order.
- Produces: Python `espressomd.shapes.VoxelMask(file, dataset, voxel_size, origin=(0,0,0), halo=2)`.
- Loader is idempotent: it stores the local box corners and periodicity of the last load and returns early when they are unchanged. Every call method of the script `VoxelMask` first calls `shape()->on_local_geometry_change()`, so a shape outside any constraint is never stale.
- Loader algorithm: for each axis, map every stored-region index to `mask_index`. Group consecutive runs of increasing mask indices into segments. For each segment triple, read one hyperslab and scatter it. Then call `set_local_data`. Geometry comes from `::System::get_system()`: `box_geo` (box_l, periodicity) and `local_geo` (`my_left()`, `my_right()`). Wrap construction in `context()->parallel_try_catch`.

- [ ] **Step 1: Write the failing tests** in `shape_voxel_mask.py`. Use `@utx.skipIfMissingFeatures(["H5MD"])` and `@utx.skipIfMissingModules("h5py")`. Write files with h5py into a `tempfile.TemporaryDirectory`. The reference is a numpy brute-force distance to solid cubes, with periodic images by box length.
  - `test_distance_matches_reference`: a seeded random porous 16³ mask (30 % solid), box 16, halo 3. For 200 seeded positions in pore, `calc_distance` matches the reference to 1e-12. For positions in solid voxels, `is_inside` is True.
  - `test_orientation`: a 8×12×16 mask with one solid voxel at index (1, 5, 11), box (8, 12, 16). `is_inside([1.5, 5.5, 11.5])` is True and `is_inside([11.5, 5.5, 1.5])` is False.
  - `test_mask_smaller_than_box`: a mask spanning x ∈ [0, 8) in a box of 12, all solid. `is_inside([10, 1, 1])` is False. `calc_distance` at x = 10 equals 2.0.
  - `test_mask_larger_than_box`: a mask of 20 voxels in x, box 16, solid only at x index 18. Every position is pore.
  - `test_local_memory`: for each node grid with product equal to `system.cell_system.get_state()["n_nodes"]`, set `system.cell_system.node_grid`. Each entry of `get_local_voxel_counts()` equals the product of stored-region extents from the Task 1 formula, computed in Python from the node grid. With more than 1 rank, each entry is smaller than the full mask voxel count. `get_local_memory_bytes()[i] <= counts[i] / 8 + 64`.
  - `test_rasterize_raises`: `rasterize(...)` raises `RuntimeError`.
  - `test_missing_file_raises`: a nonexistent path raises `RuntimeError` whose message contains the path.
  - `test_pickle`: `pickle.loads(pickle.dumps(shape))` gives identical `calc_distance` at 10 positions.

- [ ] **Step 2: Register and verify failure.**
`python_test(FILE shape_voxel_mask.py MAX_NUM_PROC 4)`, plus variants with `MAX_NUM_PROC 3 SUFFIX 3_cores`, `MAX_NUM_PROC 2 SUFFIX 2_cores` and `MAX_NUM_PROC 1 SUFFIX 1_core`.
Run: `make -C build -j8 && ctest --test-dir build -R shape_voxel_mask --output-on-failure`
Expected: FAIL with `AttributeError: module 'espressomd.shapes' has no attribute 'VoxelMask'`.

- [ ] **Step 3: Implement** the script interface class, the registration, the Python class and the docs.

- [ ] **Step 4: Run.** Same command. Expected: PASS for all four rank counts.

- [ ] **Step 5: Commit.** `git commit -m "shapes: Expose distributed VoxelMask to Python"`

---

### Task 4: Constraint coupling and lifecycle hooks

**Files:**
- Modify: `src/core/constraints/Constraint.hpp` (add `virtual void on_local_geometry_change() {}`), `src/core/constraints/ShapeBasedConstraint.hpp` (override: forward to `m_shape->on_local_geometry_change()`), `src/core/constraints/Constraints.hpp` (add `on_node_grid_change()` and `on_periodicity_change()`, both looping over constraints), `src/core/system/System.cpp` (`System::on_node_grid_change` calls `constraints->on_node_grid_change()` after `update_local_geo()`; `System::on_periodicity_change` calls `constraints->on_periodicity_change()`)
- Test: `testsuite/python/shape_voxel_mask.py` (new class `VoxelMaskConstraint`)

**Interfaces:**
- Consumes: `Shapes::Shape::on_local_geometry_change()` (Task 1), the loader set by the script interface (Task 3).

- [ ] **Step 1: Write the failing tests.** WCA between particle type 0 and constraint type 1 (σ = 1, ε = 1, voxel = 1). Skin 0.4, halo 3. Reference forces come from the numpy reference distance and the analytic WCA force.
  - `test_forces_match_reference`: the porous mask from Task 3 and 100 seeded particles in pore with reference distance > 0.6. After `integrate(0, recalc_forces=True)`, forces match the reference to 1e-10.
  - `test_forces_after_node_grid_change`: run the previous check for every node grid of the current rank count. Forces are identical across grids.
  - `test_periodicity_change`: a slab mask with solid at x index 15 only, box 16, and a particle at x = 0.5. Periodic: force (+x) is nonzero. After `system.periodicity = [False, True, True]` the force is zero. Restoring periodicity restores the force.
  - `test_force_zero_far_from_wall`: a particle more than 2 voxels from any solid. Force is exactly `[0, 0, 0]`.

- [ ] **Step 2: Run; verify failure.**
Run: `ctest --test-dir build -R shape_voxel_mask --output-on-failure`
Expected: `test_forces_after_node_grid_change` fails with a "position outside the local region" error, and `test_periodicity_change` fails because the force stays nonzero.

- [ ] **Step 3: Implement the hooks.**

- [ ] **Step 4: Run** the new tests plus `constraint_shape_based` and `shapes`.
Run: `ctest --test-dir build -R "shape_voxel_mask|constraint_shape_based|^shapes" --output-on-failure`
Expected: PASS.

- [ ] **Step 5: Commit.** `git commit -m "constraints: Reload rank-local shape data on topology changes"`

---

### Task 5: Rank-local LB boundaries from shapes

**Files:**
- Modify: `src/walberla_bridge/src/utils/boundary.hpp` (add `set_boundary_from_predicate`), `src/walberla_bridge/include/walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp`, `src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp`, `src/walberla_bridge/src/lattice_boltzmann/LBBoundaryAccess.impl.hpp`, `src/script_interface/walberla/LBFluid.cpp`, `src/python/espressomd/lb.py`, `doc/sphinx/lb.rst` (a note in the boundaries section)
- Test: `testsuite/python/shape_voxel_mask.py` (new class `VoxelMaskLB`, `@utx.skipIfMissingFeatures(["WALBERLA"])`)

**Interfaces:**
- Produces, in `walberla`:

```cpp
template <class BoundaryModel, class DataType>
void set_boundary_from_predicate(
    BoundaryModel &boundary, LatticeWalberla const &lattice,
    std::function<bool(Utils::Vector3i const &)> const &is_boundary,
    DataType const &value);
```

  Same block and ghost-layer loops as `set_boundary_from_grid`. The predicate receives the unfolded global node index, which may be negative or at least the grid size.
- Produces in `LBWalberlaBase`: `virtual void update_boundary_from_predicate(std::function<bool(Utils::Vector3i const &)> const &is_boundary, Utils::Vector3d const &velocity) = 0;`. The implementation follows `update_boundary_from_shape` step by step.
- Produces the `LBFluid` call method `"add_boundary_from_shape_local"`, with params `shape` (a `Shapes::Shape` object) and `velocity` (3 floats, lattice units). Before evaluating, call `shape->on_local_geometry_change()` (a no-op unless stale). The predicate is `shape->is_inside(box_geo.folded_position((node + 0.5) * agrid))`. After the call, run `get_system().on_lb_boundary_conditions_change()`, as `clear_boundaries` does.
- In `lb.py`, `add_boundary_from_shape`:
  - Velocity of shape `(3,)`: call `"add_boundary_from_shape_local"`.
  - Array velocity with a `VoxelMask`: raise `ValueError`.
  - Array velocity with any other shape: keep the existing path.
  - `VoxelMask`: `warnings.warn` if `abs(voxel_size - agrid) > 1e-10 * agrid`, or if `origin / agrid` is not an integer within 1e-10.

- [ ] **Step 1: Write the failing tests.**
  - `test_boundary_flags_match_mask`: the porous 16³ mask, box 16, agrid 1. `lbf[:, :, :].is_boundary` equals the mask (as bool).
  - `test_boundary_flags_all_node_grids`: the same check for every node grid of the current rank count. Create the LB fluid after setting the grid, and remove it in `tearDown`.
  - `test_alignment_warning`: origin (0.5, 0, 0) gives `assertWarns(UserWarning)`.
  - `test_array_velocity_raises`: array velocity raises `ValueError`.

- [ ] **Step 2: Run; verify failure.**
Run: `ctest --test-dir build -R shape_voxel_mask --output-on-failure`
Expected: `VoxelMaskLB` fails with `RuntimeError("VoxelMask cannot be rasterized globally")`.

- [ ] **Step 3: Implement** the bridge helper, the virtual and its implementation, the script call method, the Python routing and the docs.

- [ ] **Step 4: Run** the new tests and all LB boundary regression tests.
Run: `make -C build -j8 && ctest --test-dir build -R "shape_voxel_mask|lb_boundary|lb_shear|ek_boundary" --output-on-failure`
Expected: PASS. The existing `lb_boundary*` tests now exercise the local path for every shape.

- [ ] **Step 5: Commit.** `git commit -m "lb: Set boundaries from shapes rank-locally"`

---

### Task 6: Format and full check

- [ ] **Step 1:** Run `maintainer/format` on every changed file. Expected: no diff afterwards.
- [ ] **Step 2:** Run `ctest --test-dir build -j8 -L "python|unit" --output-on-failure` (or the full `ctest` if those labels don't exist). Expected: no new failures compared with a run on the base commit `a3a15aa4ac`.
- [ ] **Step 3:** Commit any formatting changes: `git commit -m "Format"`.
