# Third-party headers

## pystencils_runtime

Runtime headers of [pystencils](https://i10git.cs.fau.de/pycodegen/pystencils)
v2.1 (GPL-3.0-or-later), included by the kernels in `*/generated_kernels/`.
Only the headers the committed kernels include, directly or transitively,
are copied verbatim from the pystencils Python package. They must match the
pystencils version used to generate the kernels
(see `maintainer/walberla_kernels/Readme.md`).

Refresh them from the code generation virtual environment after regenerating
kernels or changing the pystencils version:

```sh
python3 $(git rev-parse --show-toplevel)/maintainer/walberla_kernels/vendor_pystencils_runtime.py
```
