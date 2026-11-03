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
import warnings

import packaging.specifiers
import packaging.version
import pystencils as ps
import pystencilssfg
import pystencilssfg.ir.entities
import pystencilssfg.ir.postprocessing
import sympy
import lbmpy
import sweepgen
import sweepgen.api
import sweepgen.build_config
import sweepgen.core.blockforest_extraction

PRECISIONS = ("single", "double")
TARGETS = ("cpu", "avx", "cuda")

# Lower bounds of the toolchain: older versions lack features that the
# generators use, so they abort.
TOOLCHAIN_REQUIREMENTS = {
    ps: ">=2.1.0",
    lbmpy: ">=2.1.1",
    pystencilssfg: ">=0.3.0",
}
# The versions the committed kernels were generated and verified with. Any
# other version of the toolchain only produces a warning.
TESTED_TOOLCHAIN = {
    ps: "2.1.0",
    lbmpy: "2.1.1",
    pystencilssfg: "0.3.0",
}
# The abbreviated commit hash of a sweepgen version, ``7.3+8.gecbe1315``
SWEEPGEN_COMMIT_RE = re.compile(r"\.g([0-9a-f]+)$")
CMAKELISTS = pathlib.Path(__file__).resolve().parents[2] / "CMakeLists.txt"
KERNEL_NAMESPACE = "walberla::pystencils"


class NameOrderedSet(set):
    """A set whose iteration follows the names of its elements."""

    def __iter__(self):
        return iter(sorted(set.__iter__(self), key=lambda x: x.name))


def _name_ordered_property(cls, name):
    """Make the property ``cls.name`` return a `NameOrderedSet`."""
    getter = getattr(cls, name).fget

    def ordered(self):
        return NameOrderedSet(getter(self))

    ordered.name_ordered = True
    setattr(cls, name, property(ordered))


def order_generated_statements():
    """
    Make the order of the statements of the generated sources independent
    of the run.

    sweepgen and pystencils-sfg emit the extraction of the fields of a
    kernel, of its scalar properties and of the components of its vector
    properties in the iteration order of sets of fields and variables. The
    hash of these objects contains the hash of their data type, which depends
    on the address of its class, so the order changes between runs, even with
    a fixed ``PYTHONHASHSEED``. The iteration of these sets is made to follow
    the names of the elements (`NameOrderedSet`), which makes the generated
    sources byte-identical between runs.

    The patch replaces four internal names, which are not a stable API:
    ``SfgKernelHandle.fields``, ``SfgKernelHandle.scalar_parameters`` and
    ``PostProcessingContext.live_variables`` of pystencils-sfg, and
    ``BlockforestParamExtraction.filter_params`` of sweepgen. A renamed
    attribute fails at import with an ``AttributeError``; a change in the
    way sweepgen or pystencils-sfg iterate the sets shows up in the
    determinism tests of the generators. The patch is applied once, also if
    the module is imported again or reloaded.
    """
    entities = pystencilssfg.ir.entities
    handle = entities.SfgKernelHandle
    if getattr(handle.fields.fget, "name_ordered", False):
        return
    _name_ordered_property(handle, "fields")
    _name_ordered_property(handle, "scalar_parameters")
    _name_ordered_property(
        pystencilssfg.ir.postprocessing.PostProcessingContext,
        "live_variables")

    extraction = (
        sweepgen.core.blockforest_extraction.BlockforestParamExtraction)
    filter_params = extraction.filter_params

    def ordered_filter_params(self, params):
        remaining, need_blockforest = filter_params(self, params)
        return NameOrderedSet(remaining), need_blockforest

    extraction.filter_params = ordered_filter_params


order_generated_statements()


def check_versions():
    """
    List the toolchain modules that are older than the lower bounds in
    ``TOOLCHAIN_REQUIREMENTS``. Warn about those that differ from the
    version in ``TESTED_TOOLCHAIN``.
    """
    errors = []
    for module, requirement in TOOLCHAIN_REQUIREMENTS.items():
        installed = packaging.version.Version(module.__version__)
        if installed not in packaging.specifiers.SpecifierSet(requirement):
            errors.append(
                f"{module.__name__} version {module.__version__} doesn't "
                f"satisfy {requirement}")
        elif installed != packaging.version.Version(
                TESTED_TOOLCHAIN[module]):
            warnings.warn(
                f"{module.__name__} version {module.__version__} is "
                f"untested (tested: {TESTED_TOOLCHAIN[module]}): run the "
                f"test suite and compare the regenerated kernels",
                stacklevel=2)
    return errors


def walberla_commit(cmake_text=None):
    """
    The waLBerla commit ESPResSo builds against: the ``GIT_TAG`` of the
    ``FetchContent_Declare(walberla ...)`` block of the top-level
    ``CMakeLists.txt`` (or of the text ``cmake_text``).
    """
    if cmake_text is None:
        cmake_text = CMAKELISTS.read_text(encoding="utf-8")
    block = re.search(r"FetchContent_Declare\(\s*walberla\b.*?^\s*\)",
                      cmake_text, re.DOTALL | re.MULTILINE)
    tag = block and re.search(r"^\s*GIT_TAG\s+([0-9a-f]{7,40})\b",
                              block.group(0), re.MULTILINE)
    if not tag:
        raise RuntimeError(
            "couldn't find the GIT_TAG commit of "
            "FetchContent_Declare(walberla ...) in CMakeLists.txt")
    return tag.group(1)


def check_sweepgen_commit(version, commit):
    """
    List the problem if the sweepgen ``version`` (``X+N.g<abbrev>``) doesn't
    come from the waLBerla commit ``commit``.
    """
    match = SWEEPGEN_COMMIT_RE.search(version)
    if match and commit.startswith(match.group(1)):
        return []
    return [f"sweepgen version {version} doesn't come from waLBerla commit "
            f"{commit[:8]}: reinstall it from "
            f"<build-dir>/_deps/walberla-src/sweepgen after configuring "
            f"ESPResSo with that commit"]


def check_toolchain():
    """
    Abort with all problems of the code generation toolchain: versions
    below the lower bounds and a sweepgen that doesn't come from the waLBerla
    commit of ESPResSo. Warn about untested versions.

    Missing internal APIs need no check of their own: the generators fail
    with an ``ImportError`` or ``AttributeError`` when they use them, and
    the test suite catches changed behaviour.
    """
    errors = check_versions()
    errors += check_sweepgen_commit(sweepgen.__version__, walberla_commit())
    if errors:
        raise RuntimeError(
            "unsupported code generation toolchain:\n  "
            + "\n  ".join(errors))


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
    needed when the toolchain in ``TESTED_TOOLCHAIN`` changes.

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
    check_toolchain()
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
