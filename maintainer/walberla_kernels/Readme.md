# Automated kernel generation with waLBerla

The scripts in this directory generate the kernels for lattice-based algorithms.
The generated kernels are committed to the repository; CMake never runs the
code generators.

## Migration status

All lattice-Boltzmann kernels are generated with
[sweepgen](https://walberla.net) and
[pystencils-sfg](https://pycodegen.pages.i10git.cs.fau.de/pystencils-sfg/).
The electrokinetics kernels have not been migrated yet: they are still the
committed output of pystencils 1.4 / `pystencils_walberla` and Jinja
templates.

| Kernel family | Generator | Toolchain |
|---|---|---|
| LB `InitialPDFsSetter` | `generate_lb_kernels.py --kernels init` | sweepgen |
| LB `StreamCollideSweepThermalized` | `generate_lb_kernels.py --kernels stream_collide_thermalized` | sweepgen |
| LB `StreamCollideSweepAthermal` | `generate_lb_kernels.py --kernels stream_collide_athermal` | sweepgen |
| LB `StreamCollideSweepLeesEdwards` | `generate_lb_kernels.py --kernels stream_collide_lees_edwards` | sweepgen |
| LB `UpdateVelFromPDF` | `generate_lb_kernels.py --kernels vel_update` | sweepgen |
| LB CPU `FieldAccessors` | `generate_lb_kernels.py --kernels field_accessors` | sweepgen |
| LB CUDA `FieldAccessors` | `generate_lb_kernels.py --kernels field_accessors_cuda` | sweepgen |
| LB `PackInfoPdf` (streaming-aware) | `generate_lb_kernels.py --kernels pack_info_pdf` | sweepgen |
| LB `PackInfoVec` | `generate_lb_kernels.py --kernels pack_info_vec` | sweepgen |
| LB `DynamicUBB` (boundary) | `generate_lb_kernels.py --kernels dynamic_ubb` | sweepgen |
| EK kernels | not regenerable at the moment | pystencils 1.4 (committed output) |

To regenerate an EK kernel, use the generator scripts of commit
`fd77c551db`, the last commit before this migration
(e.g. `git worktree add ../espresso-old-codegen fd77c551db`), together with
the pystencils 1.4 environment described in that commit's Readme.
`generate_ek_kernels.py` in the current tree no longer works with either
toolchain: it depends on helpers that were removed from `pystencils_espresso.py`
(`generate_config`, `precision_rng_modulo`) and on `pystencils_walberla`. The
modules that only it uses (`custom_additional_extensions.py`,
`walberla_ek_generation.py`, `code_generation_context.py`, `ekin.py` and the
Jinja templates in `templates/`) stay until the EK kernels are migrated.

## Code generation environment

The sweepgen kernels need pystencils, lbmpy, pystencils-sfg and the sweepgen
package of the exact waLBerla commit ESPResSo is built against
(the waLBerla `GIT_TAG` in the `FetchContent_Declare(walberla ...)` block of the top-level `CMakeLists.txt`):

```sh
python3 -m venv sweepgen
source sweepgen/bin/activate
python3 -m pip install -c "$(git rev-parse --show-toplevel)/requirements.txt" \
    pystencils lbmpy pystencilssfg sympy packaging clang-format==19.1.7
python3 -m pip install -e <build-dir>/_deps/walberla-src/sweepgen   # e.g. build_sweepgen
deactivate
```

## Generating kernels

Every kernel is written to its own file pair, named after its C++ class,
so a single kernel can be regenerated without touching the others:

```sh
source sweepgen/bin/activate
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/lattice_boltzmann/generated_kernels/
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py \
    --kernels init --precision double --target cuda
```

`--kernels`, `--precision` (`single`, `double`) and `--target` (`cpu`, `avx`,
`cuda`) all accept `all`, which is also the default. Combinations a kernel
does not support are skipped with a notice.

Format the generated files afterwards:

```sh
function format_kernels {
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cu  -style "{Language: Cpp, ColumnLimit: 0}"
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cuh -style "{Language: Cpp}"
}
format_kernels
```

The generator doesn't let pystencils-sfg format the files, even if clang-format
is on `PATH`, so the output doesn't depend on the environment.

The output is deterministic, except for the order of the field-extraction
statements in the `.cpp` and `.cu` files of the `init`, `vel_update` and
`stream_collide_*` kernels, which pystencils-sfg 0.3 does not keep stable
between runs. All headers, both `field_accessors` kernels and all files of the
`pack_info_*` and `dynamic_ubb` kernels are byte-identical for any
`PYTHONHASHSEED`.

## Architecture

* `sfg_driver.py`: toolchain version check, precision/target handling,
  one named file pair per kernel (not formatted, see above).
* `sfg_driver.generate_sweep()`: emits a sweep as a plain sweepgen class in
  `walberla::pystencils` (see "Generated sweep classes"). The classes of the
  LB sweeps keep the names of the old `pystencils_walberla` classes, but not
  their interface.
* `sfg_driver.patch_kernel_source()`: regex workarounds for pystencils 2.1.0
  code generation issues in the generated `.cpp`/`.cu` files; re-check them
  whenever the pystencils version changes.
* `lees_edwards.py`: Lees-Edwards boundary layers of the stream-collide kernel.
* `pack_infos.py`: the pack infos `PackInfoPdf` and `PackInfoVec` of
  `walberla::pystencils` on the CPU and CUDA targets (there is no AVX
  variant), with the interface of the old `pystencils_walberla` classes. A
  pack info is described by the field components that it communicates in each
  direction: the populations that cross the face, edge or corner during pull
  streaming (PDFs) or all components (vector fields). Pack, unpack and local
  copy kernels are generated for each distinct set, and a class dispatches on
  the direction. The CPU classes derive from the hand-written CRTP base
  `src/walberla_bridge/src/utils/UniformCpuFieldPackInfoBase.hpp`, which the
  generated header includes by its relative path; the CUDA classes derive from
  `walberla::sweepgen::UniformGpuFieldPackInfoBase` of the waLBerla fork.
  The code relies on sweepgen internals (`PackingKernelsContext`,
  `CaptureBufferPointer`, `GhostLayerFieldExtraction`; see the module
  docstring): re-check them whenever the waLBerla commit changes.
  pystencils linearizes a buffer with the first axis of the field slowest, so
  the symbolic fields have the axes `(z, y, x)` and the buffer is x-fastest
  like in the old pack infos. The kernels of a class live in the namespace
  `<ClassName>_kernels`, so the CPU and CUDA pack infos never define the same
  function (`tests/test_pack_info_symbols.py`).
* `boundaries.py`: the boundary conditions on the links between fluid and
  boundary cells, with the interface of the old `pystencils_walberla` classes:
  `DynamicUBB` of `walberla::lbm` on the CPU and CUDA targets (there is no AVX
  variant). The generated class derives from the hand-written CRTP base
  `walberla::IndexVectorBoundary` of
  `src/walberla_bridge/src/utils/IndexVectorBoundary.hpp`, which the
  generated header includes by its relative path. The base class collects the
  links from the flag field into the `ALL`, `INNER` and `OUTER` index vectors
  of each block, copies them to the device, runs the kernel on a vector, and
  keeps one force vector for each index vector (see the documentation of the
  header for the semantics of `run`, `inner`, `outer`, `getForce` and
  `getForceVector`). The generated class provides the link kernel and the
  structs of a link (`<ClassName>IndexInfo`: the `int32_t` coordinates `x`,
  `y`, `z` and the direction `dir`, then the data of the boundary, e.g.
  `vel_0` to `vel_2` of the UBB) and of its force (`<ClassName>ForceStruct`,
  `F_0` to `F_2`). The kernel is a sparse pystencils kernel over a field of
  these structs (`CreateKernelConfig.index_field`). The generic part
  (`generate_link_boundary`) takes the assignments of the kernel and doesn't
  depend on the boundary condition, but it is LB-shaped: it requires a force
  field, and the base class hard-codes the force struct `F_0` to `F_2`. The
  boundaries of the electrokinetics (no force) need a change of the header and
  a force-optional generator, which are deferred to sub-project 3.
  `generate_dynamic_ubb` defines the UBB with lbmpy. The assignments are
  built like `create_lattice_boltzmann_boundary_kernel` of lbmpy 2.1.1 does,
  and the code relies on lbmpy, pystencils, pystencils-sfg and sweepgen
  internals (see the module docstring): re-check them whenever the toolchain
  changes. The coordinates of the links are 32 bits wide, as in the old
  kernel, and so is all the index arithmetic of the kernel (the old one used
  64 bits for the offsets in the PDF field): the PDF field of a block with
  its ghost layers must have fewer than 2^31 elements. The generated code
  doesn't check this: the bridge does (see "Index limit of the PDF field").
  The kernels of a class live in the namespace `<ClassName>_kernels`, so the
  CPU and CUDA boundaries never define the same function
  (`tests/test_boundary_symbols.py`).
  The C++ tests of the hand-written class are
  `src/walberla_bridge/tests/IndexVectorBoundary_test{,_cuda}`.
* LB field accessors (`walberla::lbm::accessor`): one function catalogue,
  two backends.
  * `accessor_equations.py`: the equations of the accessors, derived from the
    lbmpy method.
  * `accessor_catalogue.py`: the format of the catalogue. Each accessor is
    described once: its namespace, name, kind (cell, interval, whole field,
    reduction, point, list, gather, scatter), fields, kernel and the targets
    it exists on. Per-target variants carry the parameter names and the
    return type, because the API of each target is frozen.
  * `lb_accessors.py`: the LB catalogue (kernels as pystencils assignments over
    the fields of the lattice) and the entry point
    `generate_field_accessors`. The functions that exist on one target only
    are flagged in the catalogue: `Equilibrium` and `EquilibriumDistribution`
    on the CPU, `Interpolation`, `Vector::set_from_list` and the `atomicAdd`
    kernel on CUDA. The catalogue format (`accessor_catalogue.py`) and
    the backends are generic and meant to be reused for the EK accessors with
    an EK catalogue of their own, because the backends take sizes, fields and
    equations from the catalogue.
  * `accessor_backend_cpu.py`: the CPU wrappers, in a header-only
    `FieldAccessors<Precision>Precision.h`. The kernels are pystencils
    `GenericCPU` kernels defined `inline` in the header, so the bridge needs
    no source file. They live in nested namespaces `cpu_kernels` (the CUDA
    kernels in `kernels`): both have the same signatures, and a shared name
    would let the linker bind the CPU kernels to the CUDA host stubs in builds
    with CUDA (`tests/test_accessor_symbols.py`). Single-cell values are
    exchanged through the stack, interval values through `std::vector`, and
    the reductions use `AddReductionAssignment` (a sequential sum in the
    order of the cells). The header needs the pystencils runtime headers
    like the other kernels. `sfg_driver.patch_kernel_source()` doesn't process
    headers; the accessor kernels don't need it (the RNG, `sqrt`, `fmin` and
    `fmax` are not used), which `tests/test_field_accessors.py` guards.
  * `accessor_backend_cuda.py`: the CUDA wrappers, declared in a `.cuh` and
    defined in a `.cu`. The slice operations and the particle interpolation
    are pystencils CUDA kernels. The reductions fill a device buffer that is
    summed on the host in cell order, so they are deterministic. Only
    `Interpolation::add_force` scatters with `atomicAdd`, in a hand-written
    kernel.
  * `accessor_sfg.py`: pystencils-sfg helpers and the base class of the
    backends. Overloaded free functions rely on pystencils-sfg internals
    (`AllowOverloads`, `LinearBuffer`); re-check them whenever the
    pystencils-sfg version changes.
  * `tests/data/field_accessors{,_cuda}_api_*.json` hold the API expected by
    the waLBerla bridge on each target.
* `tests/`: unit tests, run with
  `python3 -m unittest discover -s tests -t .` from this directory.

The kernels include the pystencils runtime headers (`pystencils_runtime/`).
They are not copied into ESPResSo: at configure time, CMake takes them from
the pystencils package of the Python interpreter it uses
(`cmake/FindPystencils.cmake`). Building with waLBerla therefore requires
pystencils from `requirements.txt` in that Python environment.

## Hand-written C++ headers of the bridge

Two families of generated classes derive from hand-written headers in
`src/walberla_bridge/src/utils/`. The generated headers include them by their
relative path, because the include path of the bridge targets doesn't contain
`src/`. A change that keeps the interface that the generators rely on doesn't
need a regeneration.

* `UniformCpuFieldPackInfoBase.hpp`: the CRTP base of the CPU `PackInfoPdf`
  and `PackInfoVec`, the CPU counterpart of `sweepgen`'s
  `UniformGpuFieldPackInfoBase`. It implements the interface of
  `walberla::communication::UniformPackInfo` (buffers, sizes, local copy); the
  generated class provides the pack, unpack and copy kernels per direction.
* `IndexVectorBoundary.hpp`: the CRTP base of the link boundaries (`DynamicUBB`),
  with the index vectors, force vectors and the host/device handling (see
  the description of `boundaries.py` above). Its C++ tests are
  `src/walberla_bridge/tests/IndexVectorBoundary_test{,_cuda}`.

## Index limit of the PDF field

The LB kernels that compute their offsets with `int32_t` are the
single-precision AVX stream-collide kernels (`index_dtype=int32`) and the link
kernels of the boundary conditions (`DynamicUBB`, all precisions and
targets). The PDF field of a block, including its ghost layers, must
therefore have fewer than 2^31 elements:
`Q * pad16(cells_x + 2 * ghost_layers) * (cells_y + 2 * ghost_layers) *
(cells_z + 2 * ghost_layers) < 2^31`, with `Q = 19` and `pad16(n)` the
rounding up to a multiple of 16. The rounding accounts for the aligned
allocator of the AVX kernels, which pads the x extent of the allocation (and
thus the strides of the kernels) to a multiple of the SIMD alignment (at most
64 bytes, i.e. 16 single-precision elements); it is applied for every build,
overestimating the size by less than 1%.
The generated code can't check this. The constructor of `LBWalberlaImpl`
(`src/walberla_bridge/src/lattice_boltzmann/LBWalberlaImpl.hpp`) does, with
`walberla::detail::check_pdf_field_index_range`, before it allocates any field,
and throws a `std::runtime_error` (a `RuntimeError` in Python). The remedy is
more MPI ranks or blocks, or a smaller grid; the test is in
`src/walberla_bridge/tests/LBWalberlaImpl_unit_tests.cpp`.

The stream-collide sweeps share one shadow PDF field between all blocks of a
sweep object, which is cloned from the first block that is processed. All
blocks must therefore have the same size, which `LatticeWalberla` guarantees
by creating a uniform block grid; the constructor of `LBWalberlaImpl` checks
that the PDF fields of all local blocks have the same size.

## Generated sweep classes

The classes of `InitialPDFsSetter`, `UpdateVelFromPDF`,
`StreamCollideSweepThermalized`, `StreamCollideSweepAthermal` and
`StreamCollideSweepLeesEdwards` are the classes of sweepgen, in the namespace
`walberla::pystencils`: there is no wrapper around them, and the waLBerla
bridge (`src/walberla_bridge/src/lattice_boltzmann/`) calls them directly.
Their interface:

* A constructor with all properties: first the field IDs, then the scalars,
  each in alphabetical order, e.g.
  `StreamCollideSweepThermalizedDoublePrecision(forceId, pdfsId,
  block_offset, kT, omega_bulk, omega_even, omega_odd, omega_shear, seed,
  time_step)`.
* Reference accessors of the scalars (`kT()`, `omega_shear()`, `v_s()`,
  `time_step()`, ...) that the bridge uses to update a property of an existing
  object. The field IDs are read-only.
* `operator()(IBlock *)` and `runOnCellInterval(IBlock *, CellInterval)`. The
  CUDA classes take an additional `walberla::v8::sweep::exectag::GPU`, which
  carries the stream and defaults to the default stream. The stream-collide
  sweeps swap the PDF field with a shadow field (`swapShadowBuffers`) after
  each run, so only `operator()` is meaningful for them; the shadow field is
  private to the sweep object.
* The integer properties are `int64_t`, e.g. the seed and the time step of the
  RNG. The three block offsets of the RNG are one property
  `block_offset` of type `walberla::Vector3<int64_t>`.

The classes don't know which block they run on, so the bridge sets
`block_offset` to the lower corner of the global cell bounding box of the
block before it runs the thermalized sweep on that block (the RNG of a cell
depends on its global coordinates, so the noise doesn't depend on the domain
decomposition); see `StreamCollideSweepVisitor` in `LBWalberlaImpl.hpp`.

## Thermalized LB kernels

The thermal fluctuations of the stream-collide kernel use the Philox RNG of
the pystencils 2.1 runtime. Its counter layout differs from the pystencils 1.4
kernels, so a given seed produces a different (statistically equivalent)
noise realisation than in earlier ESPResSo versions.

The thermalized and athermal stream-collide steps are two classes of their
own, `StreamCollideSweepThermalized` (with the RNG) and
`StreamCollideSweepAthermal` (without; its results agree with the thermalized
kernel at `kT = 0` up to floating-point rounding). Evaluating the RNG costs
2 to 3.3 times as much per cell, and the vectorizer can't generate a branch
on `kT`, so the choice isn't made in the kernel. The bridge chooses once in
`LBWalberlaImpl::set_collision_model(kT, seed)` (the only place where `kT` is
set): the thermalized sweep if `kT > 0`, otherwise the athermal sweep. Both
are alternatives of the `CollisionModel` variant of `LBWalberlaImpl`, next to
the Lees-Edwards sweep. Only the thermalized sweep has an RNG state (the time
step counter, `get_rng_state()` and `set_rng_state()`): for `kT = 0` the
state is `None` and setting it raises an error.
