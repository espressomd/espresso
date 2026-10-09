# Third-party headers

## pystencils_runtime

Runtime headers of [pystencils](https://i10git.cs.fau.de/pycodegen/pystencils)
v2.1 (GPL-3.0-or-later), included by the kernels in `*/generated_kernels/`.
They are copied verbatim from the pystencils Python package and must match the
pystencils version used to generate the kernels
(see `maintainer/walberla_kernels/Readme.md`).

To refresh them from the code generation virtual environment:

```sh
RT=$(python3 -m pystencils.include -s)
rm -r pystencils_runtime && cp -r "${RT}/pystencils_runtime" .
```
