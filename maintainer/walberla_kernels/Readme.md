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
| Other LB kernels | not regenerable at the moment | pystencils 1.4 (committed output) |
| EK kernels | not regenerable at the moment | pystencils 1.4 (committed output) |

To regenerate a kernel that has not been migrated yet, use the generator
scripts of commit `fd77c551db`, the last commit before this migration
(e.g. `git worktree add ../espresso-old-codegen fd77c551db`), together with
the pystencils 1.4 environment described in that commit's Readme.
`generate_ek_kernels.py` and `lees_edwards.py` in the current tree no longer
work with either toolchain: they depend on helpers that were removed from
`pystencils_espresso.py`.

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
* `tests/`: unit tests, run with
  `python3 -m unittest discover -s tests -t .` from this directory.

The pystencils runtime headers used by the kernels are vendored in
`src/walberla_bridge/src/third_party/pystencils_runtime/`; refresh them
whenever the pystencils version changes.
