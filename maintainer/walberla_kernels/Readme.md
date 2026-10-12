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
| Other LB kernels (CUDA accessors, pack infos, boundaries) | not regenerable at the moment | pystencils 1.4 (committed output) |
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
}
format_kernels
```

The output is reproducible up to the order of the field-extraction statements,
which pystencils-sfg 0.3 does not keep stable between runs.

## Architecture

* `sfg_driver.py`: toolchain version check, precision/target handling,
  one named file pair per kernel.
* `compat.py`: emits each kernel as a sweepgen class in
  `walberla::pystencils::internal` plus a wrapper class with the legacy
  `pystencils_walberla` interface used by the waLBerla bridge.
* `sfg_driver.patch_kernel_source()`: regex workarounds for pystencils 2.1.0
  code generation issues in the generated `.cpp`/`.cu` files; re-check them
  whenever the pystencils version changes.
* `lees_edwards.py`: Lees-Edwards boundary layers of the stream-collide kernel.
* `accessor_equations.py`, `field_accessors.py`: the LB field accessors
  (`walberla::lbm::accessor`); the equations come from lbmpy, the header is
  composed with pystencils-sfg. Overloaded free functions rely on
  pystencils-sfg internals (`_AllowOverloads`); re-check it whenever the
  pystencils-sfg version changes. `tests/data/field_accessors_api_*.json`
  holds the accessor API expected by the waLBerla bridge.
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
