from pystencils.astnodes import LoopOverCoordinate
from pystencils.field import fields
from pystencils.data_types import type_all_numbers
from pystencils import Assignment

from lbmpy import LBMConfig, LBStencil, LBMOptimisation, Stencil, Method, ForceModel
from lbmpy.advanced_streaming.utility import get_accessor, Timestep
from lbmpy.creationfunctions import create_lb_update_rule
from lbmpy.updatekernels import create_stream_pull_with_output_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header, generate_pack_info_for_field
from lbmpy_walberla import generate_lb_pack_info

import sympy as sp


# For advanced streaming patterns (AA, EsoTwist) the timestep is seperated into Odd and Even steps. In each of these
# steps a different streaming is executed. For more common two population methods this is not the case.
# In lbmpy we indicate this with Timestep.BOTH
streaming_pattern = "pull"
accessor = get_accessor(streaming_pattern, Timestep.BOTH)

omega = 1.0  # relaxation rate of first component
shear_velocity = 0.1  # shear velocity
shear_dir = 0  # direction of shear flow
shear_dir_normal = 1  # direction normal to shear plane, for interpolation

stencil = LBStencil(Stencil.D3Q19)
q = stencil.Q
dim = stencil.D

pdfs, pdfs_tmp = fields(
    f"pdfs({q}), pdfs_tmp({q}): double[{dim}D]", layout='fzyx')
velocity_field, force_field, density_field = fields(f"velocity({dim}), force({dim}), density(1) : double[{dim}D]",
                                                    layout='fzyx')

counters = [LoopOverCoordinate.get_loop_counter_symbol(i) for i in range(dim)]
points_up = sp.Symbol('points_up')
points_down = sp.Symbol('points_down')

u_p = sp.Piecewise((1, sp.And(type_all_numbers(counters[1] <= 0, 'int'), points_down)),
                   (-1,
                    sp.And(type_all_numbers(counters[1] >= 63,
                                            'int'),
                           points_up)),
                   (0, True)) * shear_velocity

lb_config = LBMConfig(stencil=stencil, method=Method.SRT, relaxation_rate=omega, compressible=True,
                      velocity_input=velocity_field.center_vector +
                      sp.Matrix([u_p, 0, 0]),
                      force_model=ForceModel.GUO,
                      force=force_field.center_vector, kernel_type='collide_only')
lbm_opt = LBMOptimisation(symbolic_field=pdfs)

collision = create_lb_update_rule(
    lbm_config=lb_config,
    lbm_optimisation=lbm_opt)

to_insert = [s.lhs for s in collision.subexpressions
             if collision.method.first_order_equilibrium_moment_symbols[shear_dir]
             in s.free_symbols]
for s in to_insert:
    collision = collision.new_with_inserted_subexpression(s)
ma = []
for a, c in zip(collision.main_assignments, collision.method.stencil):
    if c[2] == -1:
        b = (False, True)
    elif c[2] == 1:
        b = (True, False)
    else:
        b = (False, False)
    a = Assignment(a.lhs, a.rhs.replace(points_down, b[0]))
    a = Assignment(a.lhs, a.rhs.replace(points_up, b[1]))
    ma.append(a)
collision.main_assignments = ma

stream = create_stream_pull_with_output_kernel(collision.method, pdfs, pdfs_tmp,
                                               {'velocity': velocity_field},
                                               accessor=accessor)

init = macroscopic_values_setter(collision.method, velocity=velocity_field.center_vector,
                                 pdfs=pdfs, density=density_field.center,
                                 streaming_pattern=streaming_pattern)

stencil_typedefs = {'Stencil_T': stencil}
field_typedefs = {'PdfField_T': pdfs,
                  'VectorField_T': velocity_field,
                  'ScalarField_T': density_field}

with CodeGeneration() as ctx:
    # sweeps
    generate_sweep(ctx, 'LeesEdwards_Collision', collision)
    generate_sweep(ctx, 'LeesEdwards_Stream', stream,
                   field_swaps=[(pdfs, pdfs_tmp)])
    generate_sweep(ctx, 'LeesEdwards_Setter', init)

    # communication
    generate_lb_pack_info(
        ctx,
        'LeesEdwards_PackInfo',
        stencil,
        pdfs,
        streaming_pattern=streaming_pattern)
    generate_pack_info_for_field(
        ctx, 'LeesEdwards_PackInfo_Velocity', velocity_field)

    # Info header containing correct template definitions for stencil and field
    generate_info_header(ctx, 'LeesEdwards_InfoHeader',
                         stencil_typedefs=stencil_typedefs, field_typedefs=field_typedefs)
