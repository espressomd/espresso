# Distributed voxel-mask shape

Date: 2026-10-24

## Intent

Users want to simulate particles and LB fluid inside large porous
geometries, for example a 10k x 10k x 10k solid/pore mask from a
tomography scan. The solid phase acts as a wall for particles and as a
no-slip boundary for the LB fluid.

No single MPI rank can hold such a geometry. At 1 bit per voxel, the full
mask needs 125 GB; at 8 bytes per voxel it needs 8 TB.

## Why the existing features do not work

### `PotentialField` / `ForceField`

- Python flattens the full numpy array and passes it as one
  `std::vector<double>` (`src/python/espressomd/constraints.py`,
  `_Interpolated.__init__`).
- `GlobalContext` packs the constructor parameters on rank 0 and
  broadcasts them to all ranks.
- `FieldCoupling::Fields::Interpolated` deep-copies the whole raster into
  `m_global_field` on every rank.
- Interpolation is trilinear (`bspline_3d_accumulate<2>`). Applied to a
  0/1 mask, the force is discontinuous, the barrier is finite, and the
  wall follows the voxel staircase.

### LB boundaries from shapes

- `LBFluidWalberla.add_boundary_from_shape()` calls
  `get_shape_bitmask()`, which evaluates the shape on the full LB grid in
  Python.
- The result reaches `LBWalberlaImpl::update_boundary_from_shape()` as a
  global `std::vector<int>` (4 bytes per node, 4 TB for 10k^3 nodes).
- `get_boundary_force_from_shape()` and the EK boundaries use the same
  global-raster pattern.

## Goals

1. A new shape, `VoxelMask`, backed by a solid/pore mask that is stored
   distributed across MPI ranks.
2. Particle-wall interaction through the existing `ShapeBasedConstraint`,
   using the non-bonded interaction between the particle type and the
   constraint type.
3. LB boundaries from any shape, set up rank-locally with no global
   raster, and used for `VoxelMask`.
4. A multi-rank test whose results match a single-rank run, and a test
   that checks per-rank memory.

## Non-goals

- Load balancing. Regular domain decomposition only.
- Changing `PotentialField` or `ForceField`. They keep their global
  raster.
- EK boundaries from shapes in the rank-local path. Deferred.
  The rank-local evaluation helper should be reusable for EK later.
- `get_boundary_force_from_shape()` in the rank-local path. Deferred.
- Spatially varying slip velocity in the rank-local LB path. Only a
  uniform velocity.

## Anti-goals

- Any code path that materializes the full mask, or a full-grid LB
  raster, on one rank or in Python.
- Voxel data passing through the script interface or `GlobalContext`.

## Design decisions

### Mask input (decided)

- The mask comes from an HDF5 file. Each rank reads its own hyperslab
  collectively via parallel HDF5.
- The build already requires parallel HDF5 for
  `ESPRESSO_BUILD_WITH_HDF5` (`CMakeLists.txt`). The H5MD writer in
  `src/core/io/writer/` uses it through h5xx. `VoxelMask` is only
  available with `ESPRESSO_BUILD_WITH_HDF5`.
- Python passes file path, dataset name, voxel size, origin and halo
  width. No voxel data crosses the script interface.
- My call: the dataset is 3D, indexed `[x][y][z]` in C order (x
  slowest), matching the numpy convention of `PotentialField` and the LB
  grid accessors. The mask extent is the dataset shape times the voxel
  size. Users with `[z][y][x]` stacks transpose the file beforehand.
  The reference test uses an asymmetric geometry so that a transposed
  read fails.

### Geometry and coverage (decided)

- The user decides whether the mask covers the whole box. Origin and
  voxel size are parameters.
- Where the mask reaches a box boundary, periodicity comes from the
  System state. In a periodic direction the halo wraps by the box
  length, not by the mask length. In a non-periodic direction the halo
  ends there.
- My call: the mask geometry is therefore periodic with the box period,
  and particles and LB see the same geometry. If the mask does not span
  the box in a periodic direction, the uncovered part counts as pore and
  the wrap reaches the mask's far edge only when the mask fills the box.
  The mask may also extend beyond the box; voxels beyond the box length
  are ignored. No error in either case.
- The user configures the halo width.

### Halo validation (builder's call, suggested)

- Raise an error if the halo is smaller than the largest non-bonded
  cutoff between the constraint type and any particle type.
- Check at constraint construction and when interactions change, if a
  suitable hook exists. Otherwise check at integration start.

### Outside the mask (my call)

- Positions outside the mask's extent count as pore. No wall force.
- Users add regular `Wall` constraints at the mask edges if needed.

### Distance definition (my call)

- `calculate_dist()` returns the signed distance to the nearest solid
  voxel cube within the halo. The vector points from the nearest point on
  the solid to the query position, as for the other shapes.
- Positive in the pore, negative inside the solid. Inside the solid, the
  distance is minus the distance to the nearest pore voxel within the
  halo; this gives `penetrable` and `only_positive` a meaning.
- If no solid voxel lies within the halo, return a distance larger than
  the halo width. This is safe because the halo exceeds all cutoffs.
- Flat faces stay flat. At steps the force direction changes, but the
  force magnitude stays continuous, as with `Union`.

### Local region (my call)

- Each rank stores the voxels of its local box plus the halo.
- My call: the intended interaction range is short. With a voxel about
  one particle diameter and repulsion at contact, the cutoff is about
  1.12 voxels. A brute-force search over the halo block is acceptable.
  Large cutoffs (many voxels) are a non-goal; they need no speed
  guarantee.
- No extra per-voxel data such as a distance field.
- Storage is 1 bit per voxel in memory. The file format is 1 byte per
  voxel (HDF5 has no native bit type); nonzero means solid.
- Queries outside the local region raise an error. `rasterize()` and
  other global evaluations of `VoxelMask` raise an error.
- `VoxelMask` cannot be part of a `Union` that is evaluated globally.
  The builder decides whether to forbid it or let the outside-region
  error surface.

### Lifecycle (my call)

- When the node grid or box length changes, each rank re-reads its
  hyperslab.
- A checkpoint stores only the file path and parameters. Loading a
  checkpoint re-reads the file.

### LB boundaries (decided: in scope; mechanism my call)

- LB grid spacing equals the voxel size.
- My call: the mask Python API warns when the mask origin is not aligned
  so that each LB node center falls inside exactly one voxel. The shape
  itself does not require alignment.
- New rank-local path: each rank evaluates `is_inside()` at the LB nodes
  of its own waLBerla blocks, ghost layer included, and sets the
  boundary flags directly. Only the shape object and a uniform velocity
  cross the script interface.
- The rank-local path works for all shapes. The builder decides whether
  it replaces the current global path for uniform velocity or exists
  alongside it.
- LB and MD use the same domain decomposition, so every local LB node,
  including the ghost layer, lies within the mask's local region if the
  halo is at least one voxel.

## Constraints

- C++20, MPI. Match the style of `src/shapes/` and
  `src/script_interface/shapes/`.
- New code must run with 1, 2, 3 and 4 ranks.
- Run `maintainer/format` on changed files.

## Verification

1. **Multi-rank equivalence:** a Python test with a small mask (about
   32^3 to 64^3 voxels) with porous structure. Run on 1 and on 2/3/4
   ranks. Forces, energies and LB boundary flags match the 1-rank result.
   Particle forces match to machine precision.
2. **Reference check:** on 1 rank, `VoxelMask` forces match a
   `ShapeBasedConstraint` built from a union of `Rhomboid` voxels for the
   same simple geometry, or an analytic distance for a slab mask.
3. **Per-rank memory:** a test confirms that the mask memory on each rank
   scales with local box plus halo, not with the full mask. For example,
   expose the local voxel count and compare it to the expected value for
   the node grid.
4. **Periodicity:** halo wrap is tested for periodic and non-periodic
   directions.

## Left to the builder

- Class names, file layout, Python API parameter names.
- Bit-packing layout and the nearest-solid search algorithm.
- How the shape obtains local box geometry and periodicity from the
  System.
- How to hook re-reading into node-grid and box-length changes.
- Collective HDF5 read details, including the two-part read for a
  periodic wrap.
