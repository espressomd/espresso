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
after its class in ``walberla::pystencils``, e.g.
``InitialPDFsSetterDoublePrecisionCUDA.h``.
"""

import argparse
import dataclasses
import itertools
import pathlib
import typing

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
import lbmpy.macroscopic_value_kernels
import pystencils as ps
import sympy as sp
from pystencils.sympyextensions import DynamicType
from pystencils.sympyextensions.random import Philox

import boundaries
import lb_accessors
import lees_edwards
import pack_infos
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


def cpu_gen_config(ctx):
    """Vectorize the kernels of the AVX target (off by default in pystencils).

    In single precision, 32-bit indices are used for all AVX kernels. They
    are needed by the SIMD Philox RNG (AVX has no vectors of 8 int64) and
    are harmless for the other kernels. They limit a block to about 113M
    cells, including ghost layers. Loads and stores are unaligned
    (``assume_aligned`` is not used).
    """
    if ctx.target != "avx":
        return None
    config = ps.CreateKernelConfig()
    config.cpu.vectorize.enable = True
    config.cpu.vectorize.assume_inner_stride_one = True
    if not ctx.double_precision:
        # the SIMD Philox RNG for 8 floats needs 32-bit counters:
        # AVX has no vectors of 8 int64
        config.index_dtype = "int32"
    return config


def stream_collide_config(fields, **kwargs):
    lbm_config = lbmpy.LBMConfig(stencil=STENCIL, streaming_pattern="pull",
                                 **LBM_CONFIG_KWARGS, **kwargs)
    lbm_optimisation = lbmpy.LBMOptimisation(
        symbolic_field=fields["pdfs"],
        symbolic_temporary_field=fields["pdfs_tmp"],
        field_layout=fields["pdfs"].layout,
        cse_global=True)
    return lbm_config, lbm_optimisation


def create_thermalized_update_rule(ctx, fields, method):
    """Stream-collide update rule with Philox thermal fluctuations."""
    # the RNG parameters must have the index type of the kernel
    offsets = tuple(ps.TypedSymbol(f"block_offset_{i}", DynamicType.INDEX_TYPE)
                    for i in range(3))
    rng = Philox("fluctuation_rng", ctx.data_type,
                 seed=ps.TypedSymbol("seed", DynamicType.INDEX_TYPE),
                 offsets=offsets)
    lbm_config, lbm_optimisation = stream_collide_config(
        fields, fluctuating={"temperature": sp.Symbol("kT"), "rng": rng})
    return lbmpy.create_lb_update_rule(
        lb_method=method, lbm_config=lbm_config,
        lbm_optimisation=lbm_optimisation)


def create_athermal_update_rule(fields, method):
    """Stream-collide update rule without fluctuations."""
    lbm_config, lbm_optimisation = stream_collide_config(fields)
    return lbmpy.create_lb_update_rule(
        lb_method=method, lbm_config=lbm_config,
        lbm_optimisation=lbm_optimisation)


def create_lees_edwards_update_rule(ctx, fields, method):
    """Stream-collide update rule with Lees-Edwards boundary conditions."""
    lbm_config, lbm_optimisation = stream_collide_config(fields)
    # the collision rule is built explicitly because the Lees-Edwards shift
    # modifies it before streaming is added to the update rule
    collision_rule = lbmpy.creationfunctions.create_lb_collision_rule(
        method, lbm_config=lbm_config, lbm_optimisation=lbm_optimisation)
    collision_rule = lees_edwards.add_lees_edwards_to_collision(
        collision_rule, fields["pdfs"], shear_dir_normal=1,
        dtype=ctx.data_type)
    return lbmpy.creationfunctions.create_lb_update_rule(
        collision_rule, lbm_config=lbm_config,
        lbm_optimisation=lbm_optimisation)


def generate_init(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    sfg_driver.generate_sweep(
        sfg, name, pystencils_espresso.generate_setters(method, ctx.data_type))


def generate_stream_collide_thermalized(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    sfg_driver.generate_sweep(
        sfg, name, create_thermalized_update_rule(ctx, fields, method),
        gen_config=cpu_gen_config(ctx),
        swap_fields=[(fields["pdfs"], fields["pdfs_tmp"])])


def generate_stream_collide_athermal(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    sfg_driver.generate_sweep(
        sfg, name, create_athermal_update_rule(fields, method),
        gen_config=cpu_gen_config(ctx),
        swap_fields=[(fields["pdfs"], fields["pdfs_tmp"])])


def generate_stream_collide_lees_edwards(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    sfg_driver.generate_sweep(
        sfg, name, create_lees_edwards_update_rule(ctx, fields, method),
        gen_config=cpu_gen_config(ctx),
        swap_fields=[(fields["pdfs"], fields["pdfs_tmp"])])


def generate_vel_update(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    assignments = lbmpy.macroscopic_value_kernels.macroscopic_values_getter(
        method, None, fields["velocity"], fields["pdfs"],
        use_pre_collision_pdfs=True)
    sfg_driver.generate_sweep(sfg, name, assignments)


def generate_accessors(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    lb_accessors.generate_field_accessors(sfg, ctx, method)


def generate_dynamic_ubb(sfg, ctx, name):
    fields = pystencils_espresso.generate_fields(STENCIL, ctx.data_type)
    method = create_lb_method(fields)
    boundaries.generate_dynamic_ubb(sfg, ctx, name, method, fields["pdfs"])


def generate_pack_info_pdf(sfg, ctx, name):
    pack_infos.generate_pdf_pack_info(sfg, ctx, name, STENCIL)


def generate_pack_info_vec(sfg, ctx, name):
    pack_infos.generate_vector_pack_info(sfg, ctx, name, STENCIL)


@dataclasses.dataclass(frozen=True)
class Kernel:
    stem: str
    targets: tuple
    generate: typing.Callable
    header_only: bool = False
    header_extension: str = "h"


KERNELS = {
    "init": Kernel("InitialPDFsSetter", ("cpu", "cuda"), generate_init),
    "stream_collide_thermalized": Kernel(
        "StreamCollideSweepThermalized", ("cpu", "avx", "cuda"),
        generate_stream_collide_thermalized),
    "stream_collide_athermal": Kernel(
        "StreamCollideSweepAthermal", ("cpu", "avx", "cuda"),
        generate_stream_collide_athermal),
    "stream_collide_lees_edwards": Kernel(
        "StreamCollideSweepLeesEdwards", ("cpu", "avx", "cuda"),
        generate_stream_collide_lees_edwards),
    "vel_update": Kernel("UpdateVelFromPDF", ("cpu", "cuda"),
                         generate_vel_update),
    "dynamic_ubb": Kernel("DynamicUBB", ("cpu", "cuda"), generate_dynamic_ubb),
    "pack_info_pdf": Kernel("PackInfoPdf", ("cpu", "cuda"),
                            generate_pack_info_pdf),
    "pack_info_vec": Kernel("PackInfoVec", ("cpu", "cuda"),
                            generate_pack_info_vec),
    "field_accessors": Kernel("FieldAccessors", ("cpu",), generate_accessors,
                              header_only=True),
    "field_accessors_cuda": Kernel(
        "FieldAccessors", ("cuda",), generate_accessors,
        header_extension="cuh"),
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
                name, ctx, args.output_dir, "generate_lb_kernels.py",
                header_only=kernel.header_only,
                header_extension=kernel.header_extension) as sfg:
            kernel.generate(sfg, ctx, name)
        print(f"generated {name}")


if __name__ == "__main__":
    main()
