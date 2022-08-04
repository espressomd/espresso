# Generation of kernels for Walberla

The scripts in this directory will generate the .cpp/hpp kernels for use with 
Walberla. The generated files need to be placed in `src/walberla_bridge/generated_kernels`

WARNING: The code generation sorts the arguments alphabetically by symbol name.
If you rename something, you may have to adapt the order of arguments in calling code!

The following needs to be in the Python path:

* pystencils (https://i10git.cs.fau.de/pycodegen/pystencils)
* lbmpy (https://i10git.cs.fau.de/pycodegen/lbmpy/)
* Walberla's Python components. Here the same version should be used as is used with ESPResSo.
  One way is to use the copy in ESPResSo's build directory.

## Example

```sh
export VERSION=1.0
export DEPS="${HOME}/walberla_deps"
export PYTHONPATH="${DEPS}/${VERSION}/lbmpy:${DEPS}/${VERSION}/pystencils:${DEPS}/devel/walberla/python/"
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py --single-precision
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.{cpp,h}
mv *.{cpp,h} $(git rev-parse --show-toplevel)/src/walberla_bridge/generated_kernels/
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py --single-precision
$(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.{cpp,h}
mv ReactionKernel*.{cpp,h} $(git rev-parse --show-toplevel)/src/walberla_bridge/electrokinetics/reactions/generated_kernels/
mv *.{cpp,h} $(git rev-parse --show-toplevel)/src/walberla_bridge/electrokinetics/generated_kernels/
```
