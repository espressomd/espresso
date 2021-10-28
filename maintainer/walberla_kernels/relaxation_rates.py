# This describes the mapping between LB moments and the corresponding relaxation rates
# There are 4 relaxation rates for shear, bulk, even and odd modes,
# respectively.

from lbmpy.moments import is_bulk_moment, is_shear_moment, get_order
import sympy as sp


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
