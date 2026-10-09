#
# Copyright (C) 2020-2026 The ESPResSo project
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
Generate the LB kernels of the waLBerla bridge with sweepgen and
pystencils-sfg. Each kernel is written to its own file pair, named
after the legacy class, e.g. ``InitialPDFsSetterDoublePrecisionCUDA.h``.
"""

import argparse
import dataclasses
import itertools
import pathlib
import typing

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels

import compat
import pystencils_espresso
import relaxation_rates
import sfg_driver

STENCIL = lbmpy.LBStencil(lbmpy.Stencil.D3Q19)
LBM_CONFIG_KWARGS = dict(
    weighted=True,
    compressible=True,
    delta_equilibrium=True,
    zero_centered=True,
)


def create_lb_method(fields):
    return lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=STENCIL,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(fields["force"].center_vector),
        **LBM_CONFIG_KWARGS
    )


def generate_init(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    compat.legacy_sweep(
        sfg, ctx, name,
        pystencils_espresso.generate_setters(method, ctx.data_type))


@dataclasses.dataclass(frozen=True)
class Kernel:
    stem: str
    targets: tuple
    generate: typing.Callable


KERNELS = {
    "init": Kernel("InitialPDFsSetter", ("cpu", "cuda"), generate_init),
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the waLBerla LB kernels.")
    parser.add_argument("--kernels", nargs="+", default=["all"],
                        choices=["all", *KERNELS],
                        help="Which kernels to generate")
    parser.add_argument("--precision", nargs="+", default=["all"],
                        choices=["all", *sfg_driver.PRECISIONS],
                        help="Floating-point precision(s)")
    parser.add_argument("--target", nargs="+", default=["all"],
                        choices=["all", *sfg_driver.TARGETS],
                        help="Hardware target(s)")
    parser.add_argument("--output-dir", type=pathlib.Path,
                        default=pathlib.Path("."),
                        help="Directory to write the kernels to")
    return parser.parse_args(argv)


def expand(selection, everything):
    return (list(everything) if "all" in selection
            else list(dict.fromkeys(selection)))


def main():
    args = parse_args()
    sfg_driver.check_toolchain_versions()
    for key, precision, target in itertools.product(
            expand(args.kernels, KERNELS),
            expand(args.precision, sfg_driver.PRECISIONS),
            expand(args.target, sfg_driver.TARGETS)):
        kernel = KERNELS[key]
        if target not in kernel.targets:
            print(f"skipping {key} for target {target}: not supported")
            continue
        ctx = sfg_driver.CodegenContext(precision, target)
        name = ctx.class_name(kernel.stem)
        with sfg_driver.KernelFileGenerator(
                name, ctx, args.output_dir, "generate_lb_kernels.py") as sfg:
            kernel.generate(sfg, ctx, name)
        print(f"generated {name}")


if __name__ == "__main__":
    main()
