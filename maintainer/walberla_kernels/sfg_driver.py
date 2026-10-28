#
# Copyright (C) 2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Shared infrastructure to write one waLBerla kernel per file pair
with pystencils-sfg and sweepgen.
"""

import argparse
import copy
import dataclasses
import itertools
import pathlib
import re
import typing

import packaging.specifiers
import pystencils as ps
import pystencilssfg
import sympy
import lbmpy
import sweepgen
import sweepgen.api
import sweepgen.build_config

PRECISIONS = ("single", "double")
TARGETS = ("cpu", "avx", "cuda")

TOOLCHAIN_REQUIREMENTS = {
    ps: "==2.1.0",
    lbmpy: "==2.1.1",
    pystencilssfg: "==0.3.0",
}
WALBERLA_COMMIT = "ecbe1315"
KERNEL_NAMESPACE = "walberla::pystencils"


def check_toolchain_versions():
    """Abort if the code generation toolchain is not the pinned one."""
    for module, requirement in TOOLCHAIN_REQUIREMENTS.items():
        spec = packaging.specifiers.SpecifierSet(requirement)
        assert spec.contains(module.__version__), \
            f"{module.__name__} version {module.__version__} " \
            f"doesn't match requirement {requirement}"
    assert f"g{WALBERLA_COMMIT}" in sweepgen.__version__, \
        f"sweepgen version {sweepgen.__version__} doesn't come from " \
        f"waLBerla commit {WALBERLA_COMMIT}"


@dataclasses.dataclass(frozen=True)
class CodegenContext:
    """Precision and hardware target of one generated file pair."""
    precision: str
    target: str

    def __post_init__(self):
        assert self.precision in PRECISIONS, self.precision
        assert self.target in TARGETS, self.target

    @property
    def double_precision(self):
        return self.precision == "double"

    @property
    def data_type(self):
        return ps.create_type(
            "float64" if self.double_precision else "float32")

    @property
    def is_gpu(self):
        return self.target == "cuda"

    @property
    def precision_prefix(self):
        return ("DoublePrecision" if self.double_precision
                else "SinglePrecision")

    @property
    def target_suffix(self):
        return {"cpu": "", "avx": "AVX", "cuda": "CUDA"}[self.target]

    def class_name(self, stem):
        """
        Class and file name, e.g. ``InitialPDFsSetterDoublePrecisionCUDA``.
        """
        return f"{stem}{self.precision_prefix}{self.target_suffix}"

    def build_config(self):
        """waLBerla build configuration seen by sweepgen."""
        config = sweepgen.build_config.WalberlaBuildConfig(
            c_compiler_id="GNU",
            cxx_compiler_id="GNU",
            use_double_precision=self.double_precision,
            optimize_for_localhost=False,
            mpi_enabled=True,
            openmp_enabled=not self.is_gpu,
            hip_enabled=False,
            cuda_enabled=self.is_gpu,
            likwid_enabled=False,
        )
        if self.target == "avx":
            config.target = ps.Target.X86_AVX
        return config


def index_list_field():
    """
    The index field of a sparse sweep: the cells of the sweep's
    ``walberla::v8::sweep::SparseIndexList``, as an indexed field of
    ``walberla::v8::sweep::CellIdx`` with the members ``x``, ``y`` and ``z``.
    A sparse sweep accesses its other fields with absolute coordinates,
    ``field.absolute_access(tuple(index_list[0](c) for c in "xyz"), ())``.
    """
    return ps.Field.create_generic(
        "indexList", 1, sweepgen.api.CellIdx, field_type=ps.FieldType.INDEXED)


def generate_sweep(sfg, name, assignments, gen_config=None, swap_fields=(),
                   sparse=False):
    """
    Generate the sweepgen class ``walberla::pystencils::<name>`` of a sweep
    over the cells of the fields of ``assignments``.

    ``gen_config`` is ``None`` or a ``pystencils.CreateKernelConfig``, e.g.
    the AVX settings of ``generate_lb_kernels.cpu_gen_config``.

    ``swap_fields`` is a list of ``(field, shadow_field)`` pairs that the
    sweep swaps after each run (``Sweep.swap_fields()``). The sweep keeps
    the shadow field to itself, so all blocks that it runs on must have
    fields of the same size.

    With ``sparse``, the sweep runs on the cells of a ``SparseIndexList``
    that the constructor takes as the argument ``indexList``, instead of all
    cells of the block. The assignments address the cells through
    `index_list_field()`. pystencils can't prove that the cells of the list
    are distinct and that the sweep doesn't read what it writes in another
    cell, so the independence check is skipped (``skip_independence_check``):
    the list must not contain a cell twice. A sparse sweep has no
    ``runOnCellInterval``.

    :raises ValueError: if ``sparse`` is set and ``gen_config`` already has
        an index field other than `index_list_field()` or sets
        ``skip_independence_check=False``; or if ``gen_config`` makes the
        sweep sparse (it sets an index field) and ``sparse`` is not set, or
        the other way round
    """
    if sparse:
        gen_config = (ps.CreateKernelConfig() if gen_config is None
                      else copy.copy(gen_config))
        if gen_config.index_field is not None \
                and gen_config.index_field != index_list_field():
            raise ValueError(
                f"sweep {name}: the sparse sweep uses index_list_field(), "
                f"but gen_config sets another index field")
        if gen_config.skip_independence_check is False:
            raise ValueError(
                f"sweep {name}: the sparse sweep skips the independence "
                f"check, but gen_config sets skip_independence_check=False")
        gen_config.index_field = index_list_field()
        gen_config.skip_independence_check = True
    sweep = sweepgen.Sweep(name, assignments, gen_config)
    if sweep.sparse != sparse:
        raise ValueError(
            f"sweep {name} is {'' if sweep.sparse else 'not '}sparse, but "
            f"the argument sparse is {sparse}")
    for field, shadow_field in swap_fields:
        sweep.swap_fields(field, shadow_field)
    with sfg.namespace(KERNEL_NAMESPACE):
        sfg.generate(sweep)


def toolchain_metadata(generator_script):
    return (
        f"Generated by maintainer/walberla_kernels/{generator_script}; "
        f"do not edit by hand.\n"
        f"kernel generated with pystencils v{ps.__version__}, "
        f"lbmpy v{lbmpy.__version__}, "
        f"pystencils-sfg v{pystencilssfg.__version__}, "
        f"sweepgen v{sweepgen.__version__}, "
        f"sympy v{sympy.__version__}"
    )


class KernelFileGenerator(pystencilssfg.SourceFileGenerator):
    """
    Source file generator that writes ``<basename>.<header_extension>``
    (``.h`` by default) and ``<basename>.cpp`` (``.cu`` for CUDA) instead of
    naming the files after the generator script. With ``header_only``, only
    the header is written and all functions are generated ``inline``.
    pystencils-sfg doesn't format the files, even if it finds clang-format.
    """

    def __init__(self, basename, ctx, output_dir, generator_script,
                 header_only=False, header_extension="h"):
        self._basename = basename
        self._generator_script = generator_script
        sweepgen.build_config.DEBUG.BUILD_CONFIG = ctx.build_config()
        config = pystencilssfg.SfgConfig()
        config.extensions.header = header_extension
        config.extensions.impl = "cu" if ctx.is_gpu else "cpp"
        config.header_only = header_only
        # the output must not depend on the environment: the files are
        # formatted with ``maintainer/format/clang-format.sh``
        config.clang_format.skip = True
        config.output_directory = pathlib.Path(output_dir)
        super().__init__(config, keep_unknown_argv=True)

    def _scriptname(self):
        return f"{self._basename}.py"

    def __enter__(self):
        sfg = super().__enter__()
        sfg.prelude(toolchain_metadata(self._generator_script))
        return sfg

    def __exit__(self, exc_type, exc_value, traceback):
        result = super().__exit__(exc_type, exc_value, traceback)
        if exc_type is None:
            for path in sorted(self._output_dir.glob(f"{self._basename}.*")):
                if path.suffix in (".cpp", ".cu"):
                    source = path.read_text(encoding="utf-8")
                    patched = patch_kernel_source(source)
                    if patched != source:
                        path.write_text(patched, encoding="utf-8")
        return result


def patch_kernel_source(source):
    """
    Work around code generation problems of pystencils 2.1.0 in the
    generated ``.cpp`` and ``.cu`` files. Re-check whether they are still
    needed when the toolchain in ``TOOLCHAIN_REQUIREMENTS`` changes.

    * The kernels live in a namespace nested in ``walberla::pystencils``,
      where the unqualified ``pystencils::runtime::`` of the RNG calls
      resolves to ``walberla::pystencils::runtime`` instead of the
      ``::pystencils::runtime`` of the pystencils runtime headers.
    * The vectorized RNG functions deduce their key type from both the
      ``int64_t`` seed (``long``) and the key literal (``long long``),
      which doesn't compile on LP64 platforms: pass the type explicitly.
    * Single-precision CPU kernels call the global ``sqrt()``, ``fmin()``,
      ``fmax()``, ``pow()`` and ``fabs()``, i.e. the double overloads: use
      the ``std::`` ones.
    """
    source = re.sub(r"(?<![:\w])pystencils::runtime::",
                    "::pystencils::runtime::", source)
    source = re.sub(r"(?<![\w:.>])(sqrt|fmin|fmax|pow|fabs)\(",
                    r"std::\1(", source)
    source = re.sub(r"(::x86::random::philox_float(?:32|64))\(",
                    r"\1<int64_t>(", source)
    return source


@dataclasses.dataclass(frozen=True)
class Kernel:
    """
    A kernel of a generator script: the stem of its class name, the targets
    it exists on and the function ``generate(sfg, ctx, name)`` that writes it.
    """
    stem: str
    targets: tuple
    generate: typing.Callable
    header_only: bool = False
    header_extension: str = "h"


def parse_args(description, kernels, argv=None):
    """Parse the command line of a generator script (see `generate_main`)."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--kernels", nargs="+", default=["all"],
                        choices=["all", *kernels],
                        help="Which kernels to generate")
    parser.add_argument("--precision", nargs="+", default=["all"],
                        choices=["all", *PRECISIONS],
                        help="Floating-point precision(s)")
    parser.add_argument("--target", nargs="+", default=["all"],
                        choices=["all", *TARGETS],
                        help="Hardware target(s)")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("."),
                        help="Directory to write the kernels to")
    return parser.parse_args(argv)


def expand(selection, everything):
    """Resolve ``all`` in a command line selection; drop duplicates."""
    return (list(everything) if "all" in selection
            else list(dict.fromkeys(selection)))


def generate_main(description, generator_script, kernels, argv=None):
    """
    Command line of a generator script: write each selected combination of
    kernel (a key of the dict ``kernels`` of `Kernel`), precision and target
    to its own file pair in the output directory.
    """
    args = parse_args(description, kernels, argv)
    check_toolchain_versions()
    for key, precision, target in itertools.product(
            expand(args.kernels, kernels), expand(args.precision, PRECISIONS),
            expand(args.target, TARGETS)):
        kernel = kernels[key]
        if target not in kernel.targets:
            print(f"skipping {key} for target {target}: not supported")
            continue
        ctx = CodegenContext(precision, target)
        name = ctx.class_name(kernel.stem)
        with KernelFileGenerator(
                name, ctx, args.output_dir, generator_script,
                header_only=kernel.header_only,
                header_extension=kernel.header_extension) as sfg:
            kernel.generate(sfg, ctx, name)
        print(f"generated {name}")
