# Generation of kernels for Walberla

The scripts in this directory generate the .cpp/.hpp kernels for lattice-based algorithms.

WARNING: The code generation sorts the arguments alphabetically by symbol name.
If you rename something, you may have to adapt the order of arguments in the
calling code!

The following dependencies need to be in the Python path:

* pystencils (https://i10git.cs.fau.de/pycodegen/pystencils)
* lbmpy (https://i10git.cs.fau.de/pycodegen/lbmpy/)
* waLBerla's Python components. Here the same version should be used as
  the one used to build ESPResSo. One way is to use the copy fetched in
  ESPResSo's `build/_deps/walberla-src/python/` directory.

## Example

```sh
export VERSION=1.1.1
export DEPS="${HOME}/walberla_deps"
export PYTHONPATH="${DEPS}/${VERSION}/lbmpy:${DEPS}/${VERSION}/pystencils:${DEPS}/devel/walberla/python/"
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/lattice_boltzmann/generated_kernels/
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py --single-precision
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/generated_kernels/
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py --single-precision
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
mv ReactionKernel*.{cpp,h} $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/reactions/generated_kernels/
```
