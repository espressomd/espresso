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
python3 -m venv codegen
source codegen/bin/activate
python3 -m pip install -c "$(git rev-parse --show-toplevel)/requirements.txt" \
    numpy cython sympy islpy jinja2 setuptools packaging
python3 -m pip install \
  git+https://i10git.cs.fau.de/pycodegen/lbmpy.git@d3f62364cf512cd9d93a93d89170cfea1bb2b0a1 \
  git+https://i10git.cs.fau.de/jngrad/pystencils.git@dfd203a3f6318f5fb4cfa09f9f8c12bb45463bd7
deactivate
```

The kernels can be regenerated with this shell script:

```sh
# adapt these paths to the build environment
source codegen/bin/activate
export PYTHONPATH="$(realpath build/_deps/walberla-src/python/)"

# convenience functions
function generate_lb_kernels {
  python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_lb_kernels.py $@
}
function generate_ek_kernels {
  python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/generate_ek_kernels.py $@
}
function format_kernels {
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.h
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cpp -style "{Language: Cpp, ColumnLimit: 0}"
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cu  -style "{Language: Cpp, ColumnLimit: 0}"
  $(git rev-parse --show-toplevel)/maintainer/format/clang-format.sh -i *.cuh -style "{Language: Cpp}"
}

# LB kernels
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/lattice_boltzmann/generated_kernels/
generate_lb_kernels
generate_lb_kernels --single-precision
generate_lb_kernels --gpu
generate_lb_kernels --gpu --single-precision
format_kernels

# EK kernels
cd $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/generated_kernels/
generate_ek_kernels
generate_ek_kernels --single-precision
generate_ek_kernels --gpu
generate_ek_kernels --gpu --single-precision
format_kernels
mv ReactionKernel*.{cpp,h,cu} $(git rev-parse --show-toplevel)/src/walberla_bridge/src/electrokinetics/reactions/generated_kernels/
```

The code generation is not deterministic, therefore the list of changes might
be quite large. If you only adapted a few lines in a specific template file,
then you only need to commit the corresponding output files.

WARNING: The code generation sorts the arguments alphabetically by symbol name.
If you rename something, you may have to adapt the order of arguments in the
calling code!
