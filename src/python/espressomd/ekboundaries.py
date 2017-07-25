from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class EKBoundaries(ScriptInterfaceHelper):
    _so_name = "LBBoundaries::LBBoundaries"

    def add(self, *args, **kwargs):
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
        self.call_method("remove", ekboundary=ekboundary)


@script_interface_register
class EKBoundary(ScriptInterfaceHelper):
    _so_name = "LBBoundaries::LBBoundary"
