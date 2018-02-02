from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class EKBoundaries(ScriptInterfaceHelper):
    """
    Creates a set of electrokinetics boundaries.
    """
    _so_name = "LBBoundaries::LBBoundaries"

    def add(self, *args, **kwargs):
        """
        Adds a boundary to the set.
        Either a valid boundary is an argument,
        or a valid set of parameters to create a boundary.
        """
        if len(args) == 1:
            if isinstance(args[0], EKBoundary):
                ekboundary = args[0]
            else:
                raise TypeError(
                    "Either a EKBoundary object or key-value pairs for the parameters of a EKBoundary object need to be passed.")
        else:
            ekboundary = EKBoundary(**kwargs)
        self.call_method("add", object=ekboundary)
        return ekboundary

    def remove(self, ekboundary):
        """
        Removes a boundary from the set.
        
        Parameters
        ----------
        ekboundary : :obj: `EKBoundary`
                     The boundary to be removed from the set.
        """
        self.call_method("remove", ekboundary=ekboundary)


@script_interface_register
class EKBoundary(ScriptInterfaceHelper):
    """
    Creates a EK boundary.
    """
    _so_name = "LBBoundaries::LBBoundary"
