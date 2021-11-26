# This describes the mapping between LB moments and the corresponding relaxation rates
# There are 4 relaxation rates for shear, bulk, even and odd modes,
# respectively.

# Original source:
# https://i10git.cs.fau.de/pycodegen/lbmpy/-/blob/0e7962be84613466e6842f37111c571db8183b3d/lbmpy_tests/test_fluctuating_lb.py#L25-47

from lbmpy.moments import is_bulk_moment, is_shear_moment, get_order
import sympy as sp


def rr_getter(moment_group):
    """Maps a group of moments to a relaxation rate (shear, bulk, even, odd)
    in the 4 relaxation time thermalized LB model or 0 for conserved modes.
    """
    is_shear = [is_shear_moment(m, 3) for m in moment_group]
    is_bulk = [is_bulk_moment(m, 3) for m in moment_group]
    order = [get_order(m) for m in moment_group]
    assert min(order) == max(order)
    order = order[0]

    if order < 2:
        return [0] * len(moment_group)
    elif any(is_bulk):
        assert all(is_bulk)
        return [sp.Symbol("omega_bulk")] * len(moment_group)
    elif any(is_shear):
        assert all(is_shear)
        return [sp.Symbol("omega_shear")] * len(moment_group)
    elif order % 2 == 0:
        assert order > 2
        return [sp.Symbol("omega_even")] * len(moment_group)
    else:
        return [sp.Symbol("omega_odd")] * len(moment_group)
