from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class ComFixed(ScriptInterfaceHelper):
    """Fix the center of mass of specific types.

    Subtracts mass-weighted fraction of the total
    force action on all particles of the type from
    the particles after each force calculation. This
    keeps the center of mass of the type fixed iff
    the total momentum of the type is zero.

    Parameters
    ----------
    types : array_like
        List of types of which the center of mass
        should be fixed.
    """

    _so_name = "ComFixed"
    _so_creation_policy = "GLOBAL"
