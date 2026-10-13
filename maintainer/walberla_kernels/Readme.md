# Automated kernel generation with waLBerla

The scripts in this directory generate the kernels for lattice-based algorithms.
The generated kernels are committed to the repository; CMake never runs the
code generators.

## Migration status

The kernels are being migrated from pystencils 1.4 / `pystencils_walberla`
and Jinja templates to [sweepgen](https://walberla.net) and
[pystencils-sfg](https://pycodegen.pages.i10git.cs.fau.de/pystencils-sfg/).

| Kernel family | Generator | Toolchain |
|---|---|---|
| LB `InitialPDFsSetter` | `generate_lb_kernels.py --kernels init` | sweepgen |
| LB `StreamCollideSweepThermalized` | `generate_lb_kernels.py --kernels stream_collide_thermalized` | sweepgen |
| LB `StreamCollideSweepLeesEdwards` | `generate_lb_kernels.py --kernels stream_collide_lees_edwards` | sweepgen |
| LB `UpdateVelFromPDF` | `generate_lb_kernels.py --kernels vel_update` | sweepgen |
| LB CPU `FieldAccessors` | `generate_lb_kernels.py --kernels field_accessors` | sweepgen |
| LB CUDA `FieldAccessors` | `generate_lb_kernels.py --kernels field_accessors_cuda` | sweepgen |
| Other LB kernels (pack infos, boundaries) | not regenerable at the moment | pystencils 1.4 (committed output) |
| EK kernels | not regenerable at the moment | pystencils 1.4 (committed output) |

To regenerate a kernel that has not been migrated yet, use the generator
scripts of commit `fd77c551db`, the last commit before this migration
(e.g. `git worktree add ../espresso-old-codegen fd77c551db`), together with
the pystencils 1.4 environment described in that commit's Readme.
`generate_ek_kernels.py` in the current tree no longer works with either
toolchain: it depends on helpers that were removed from `pystencils_espresso.py`.

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
between runs. All headers and both `field_accessors` kernels are byte-identical
for any `PYTHONHASHSEED`.

## Architecture

* `sfg_driver.py`: toolchain version check, precision/target handling,
  one named file pair per kernel (not formatted, see above).
* `compat.py`: emits each kernel as a sweepgen class in
  `walberla::pystencils::internal` plus a wrapper class with the legacy
  `pystencils_walberla` interface used by the waLBerla bridge.
* `sfg_driver.patch_kernel_source()`: regex workarounds for pystencils 2.1.0
  code generation issues in the generated `.cpp`/`.cu` files; re-check them
  whenever the pystencils version changes.
* `lees_edwards.py`: Lees-Edwards boundary layers of the stream-collide kernel.
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
    no source file. Single-cell values are exchanged through the stack,
    interval values through `std::vector`, and the reductions use
    `AddReductionAssignment` (a sequential sum in the order of the cells).
    The header needs the pystencils runtime headers like the
    other kernels. `sfg_driver.patch_kernel_source()` doesn't process
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

## Thermalized LB kernels

The thermal fluctuations of the stream-collide kernel use the Philox RNG of
the pystencils 2.1 runtime. Its counter layout differs from the pystencils 1.4
kernels, so a given seed produces a different (statistically equivalent)
noise realisation than in earlier ESPResSo versions. When `kT` is zero, the
kernel wrapper runs a second, generated kernel without RNG; its results agree
with the thermalized kernel at `kT = 0` up to floating-point rounding.
