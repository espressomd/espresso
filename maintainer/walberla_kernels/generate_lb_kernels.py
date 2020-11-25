import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule, create_mrt_orthogonal, force_model_from_string
from lbmpy.moments import is_bulk_moment, is_shear_moment, get_order
from lbmpy.stencils import get_stencil
from pystencils_walberla import CodeGeneration
from lbmpy_walberla import generate_lattice_model

with CodeGeneration() as ctx:
    kT = sp.symbols("kT")
    force_field = ps.fields("force(3): [3D]", layout='zyxf')

    def rr_getter(moment_group):
        """Map groups of LB moments to relaxation rate symbols 
        (omega_shear/bulk/odd/even) or 0 for conserved modes.
        """
        is_shear = [is_shear_moment(m, 3) for m in moment_group]
        is_bulk = [is_bulk_moment(m, 3) for m in moment_group]
        order = [get_order(m) for m in moment_group]
        assert min(order) == max(order)
        order = order[0]

        if order < 2:
            return 0
        elif any(is_bulk):
            assert all(is_bulk)
            return sp.Symbol("omega_bulk")
        elif any(is_shear):
            assert all(is_shear)
            return sp.Symbol("omega_shear")
        elif order % 2 == 0:
            assert order > 2
            return sp.Symbol("omega_even")
        else:
            return sp.Symbol("omega_odd")

    method = create_mrt_orthogonal(
        stencil=get_stencil('D3Q19'),
        compressible=True,
        weighted=True,
        relaxation_rate_getter=rr_getter,
        force_model=force_model_from_string(
            'schiller', force_field.center_vector)
    )
    # generate unthermalized LB
    collision_rule_unthermalized = create_lb_collision_rule(
        method,
        optimization={'cse_global': True}
    )
    generate_lattice_model(
        ctx,
        'MRTLatticeModel',
        collision_rule_unthermalized)

    # generate thermalized LB
    collision_rule_thermalized = create_lb_collision_rule(
        method,
        fluctuating={
            'temperature': kT,
            'block_offsets': 'walberla',
        },
        optimization={'cse_global': True}
    )
    generate_lattice_model(
        ctx,
        'FluctuatingMRTLatticeModel',
        collision_rule_thermalized)