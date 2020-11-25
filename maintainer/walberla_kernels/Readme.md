# Generation of kernels for Walberla
i
The scripts in this directory will generate the .cpp/hpp kernels for use with 
Walberla. The generated files need to be placed in `src/walberla_bridge`

WARNING: The code generation sorts the arguments alphabetically by symbol name.
If you rename something, you may have to adapt the order of arguments in calling code!


Usage:
The following needs to be in the Python path:

* pystencils (https://i10git.cs.fau.de/pycodegen/pystencils)
* lbmpy (https://i10git.cs.fau.de/pycodegen/lbmpy/)
* Walberla's Python componentent. Here the same version should be used as is used with Espresso.
One way is to use the copy in Espresso's build directory.

## Example:
```
export PYTHONPATH=$HOME/pystencils:$HOME/lbmpy:$HOME/es/build-walberla/_deps/walberla-src/python/
python3 generate_lb_kernels.py
```
