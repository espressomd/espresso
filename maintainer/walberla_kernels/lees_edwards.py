
from pystencils.astnodes import LoopOverCoordinate
from pystencils.data_types import type_all_numbers
from pystencils import Assignment

from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

import sympy as sp


def velocity_offset_eqs(method, pdfs, shear_dir_normal, stencil):
    """Calculates the difference between quilibrium pdf distributions
    with (rho, u) and (rho, u+v) and applies them to out-flowing
    populations in the boundary layer. Returns an AssignmentCollection
    with one Assignment per stencil direction.
    """
    dim = len(stencil[0])

    # Placeholders indicating a population flows up or down.
    # Will be replaced later using the component of the stencil direction
    # along the shear_dir_normal.
    points_up = sp.Symbol('points_up')
    points_down = sp.Symbol('points_down')

    # Symbol for the coordinate index within the field,
    # used to identify boundary layers
    counters = [LoopOverCoordinate.get_loop_counter_symbol(
        i) for i in range(dim)]

    grid_size = sp.Symbol("grid_size", dtype=int)  # in shear_normal_dir

    # +,-1 for upper/lower boundary layers, 0 otherwise.
    # Based on symbolic counters defined above. Only becomes
    # non-zero if the corresponding points_up/down flags
    # are engaged (which is only done for out-flowing populations)
    layer_prefactor = sp.Piecewise((-1, sp.And(type_all_numbers(counters[1] <= 0, 'int'), points_down)),
                                   (1,
                                    sp.And(type_all_numbers(counters[1] >= grid_size - 1,
                                                            'int'),
                                           points_up)),
                                   (0, True))

    # Start with an equilibrium distribution for a given density and velocity
    delta_pdf_eqs = macroscopic_values_setter(
        method, sp.Symbol("dens"), [
            sp.Symbol("v_0"), sp.Symbol("v_1"), sp.Symbol("v_2")], pdfs)

    # Replace the assignments of (rho,u) by (rho, u+v) - (rho,u)
    ma = []
    for a, c in zip(delta_pdf_eqs.main_assignments, method.stencil):
        # Determine direction of the stencil component in the
        # shear_dir_normal
        if c[shear_dir_normal] == 1:
            up = True
            down = False
        elif c[shear_dir_normal] == -1:
            up = False
            down = True
        else:
            up = False
            down = False

        # Replace (rho,u) by (rho,u+v) in boundary layers
        rhs = sp.simplify(
            a.rhs -
            a.rhs.replace(
                sp.Symbol("u_0"),
                sp.Symbol("u_0") +
                layer_prefactor *
                sp.Symbol("v_s")))

        # Only engage if the population is outflowing. See layer_prefactor
        rhs = rhs.replace(points_up, up)
        rhs = rhs.replace(points_down, down)
        new_a = Assignment(a.lhs, rhs)

        ma.append(new_a)
        print(c, ma[-1])
    # Plug in modified assignments
    delta_pdf_eqs.main_assignments = ma
    return delta_pdf_eqs.main_assignments


def add_lees_edwards_to_collision(collision, pdfs, stencil, shear_dir_normal):
    # Get population shift for outflowing populations at the boundaries
    offset = velocity_offset_eqs(
        collision.method,
        pdfs,
        shear_dir_normal,
        stencil)

    ma = []
    for i, a in enumerate(collision.main_assignments):
        # Add Lees-Edwards-shift to collision main assignments
        new_a = Assignment(a.lhs, a.rhs + offset[i].rhs)
        ma.append(new_a)
    collision.main_assignments = ma
    return collision
