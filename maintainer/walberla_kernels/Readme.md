# Automated kernel generation with waLBerla

The scripts in this directory generate the kernels for lattice-based algorithms.

The following dependencies need to be in the Python path:

* pystencils (https://i10git.cs.fau.de/pycodegen/pystencils)
* lbmpy (https://i10git.cs.fau.de/pycodegen/lbmpy/)
* waLBerla's Python components. Here the same version should be used as
  the one used to build ESPResSo. One way is to use the copy fetched in
  ESPResSo's `build/_deps/walberla-src/python/` directory.

The Python dependencies can be pip installed locally with the following command:

```sh
python3 -m pip install --user -c requirements.txt numpy sympy lbmpy pystencils islpy
```

The kernels can be regenerated with this shell script:

```sh
# adapt these paths to the build environment
export VERSION=1.2
export DEPS="${HOME}/walberla_deps"
export PYTHONPATH="${DEPS}/${VERSION}/lbmpy:${DEPS}/${VERSION}/pystencils:${DEPS}/devel/walberla/python/"

# convenience functions
function generate_lb_kernels {
  python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py $@
}
function generate_ek_kernels {
  python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py $@
}
function format_lb_kernels {
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
}
function format_ek_kernels {
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
}

# LB kernels
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/lattice_boltzmann/generated_kernels/
generate_lb_kernels
generate_lb_kernels --single-precision
format_lb_kernels

# EK kernels
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/generated_kernels/
generate_ek_kernels
generate_ek_kernels --single-precision
format_ek_kernels
mv ReactionKernel*.{cpp,h} $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/reactions/generated_kernels/
```

WARNING: The code generation sorts the arguments alphabetically by symbol name.
If you rename something, you may have to adapt the order of arguments in the
calling code!
