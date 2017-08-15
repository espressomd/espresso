from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
import espressomd.code_info


if any(i in espressomd.code_info.features() for i in ["LB_BOUNDARIES", "LB_BOUNDARIES_GPU"]):
    @script_interface_register
    class LBBoundaries(ScriptInterfaceHelper):
        _so_name = "LBBoundaries::LBBoundaries"

        def add(self, *args, **kwargs):
            if len(args) == 1:
                if isinstance(args[0], LBBoundary):
                    lbboundary = args[0]
                else:
                    raise TypeError(
                        "Either a LBBoundary object or key-value pairs for the parameters of a LBBoundary object need to be passed.")
            else:
                lbboundary = LBBoundary(**kwargs)
            self.call_method("add", object=lbboundary)
            return lbboundary

        def remove(self, lbboundary):
            self.call_method("remove", lbboundary=lbboundary)


    @script_interface_register
    class LBBoundary(ScriptInterfaceHelper):
        _so_name = "LBBoundaries::LBBoundary"
        _so_bind_methods = ("get_force",)


    @script_interface_register
    class LBMovingBoundaries(ScriptInterfaceHelper):
        _so_name = "LBBoundaries::LBMovingBoundaries"

        def add(self, *args, **kwargs):
            if len(args) == 1:
                if isinstance(args[0], LBMovingBoundary):
                    lbboundary = args[0]
                else:
                    raise TypeError(
                        "Either a LBMovingBoundary object or key-value pairs for the\
                         parameters of a LBMovingBoundary object need to be passed.")
            else:
                lbboundary = LBMovingBoundary(**kwargs)
            self.call_method("add", object=lbboundary)
            return lbboundary

        def remove(self, lbboundary):
            self.call_method("remove", lbboundary=lbboundary)


    @script_interface_register
    class LBMovingBoundary(ScriptInterfaceHelper):
        _so_name = "LBBoundaries::LBMovingBoundary"
