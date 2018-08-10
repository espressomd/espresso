from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
import espressomd.code_info


if any(i in espressomd.code_info.features() for i in ["LB_BOUNDARIES", "LB_BOUNDARIES_GPU"]):
    @script_interface_register
    class LBBoundaries(ScriptInterfaceHelper):
        """
        Creates a set of lattice Boltzmann boundaries.

        """

        _so_name = "LBBoundaries::LBBoundaries"

        def __getitem__(self, key):
            return self.call_method("get_elements")[key]

        def __iter__(self):
            elements = self.call_method("get_elements")
            for e in elements:
                yield e

        def add(self, *args, **kwargs):
            """
            Adds a boundary to the set.
            Either a valid boundary is an argument,
            or a valid set of parameters to create a boundary.

            """

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
            """
            Removes a boundary from the set.

            Parameters
            ----------
            lbboundary : :obj:`LBBoundary`
                         The boundary to be removed from the set.

            """

            self.call_method("remove", object=lbboundary)

        def clear(self):
            """
            Removes all boundaries.

            """

            self.call_method("clear")

        def size(self):
            return self.call_method("size")

        def empty(self):

            return self.call_method("empty")

    @script_interface_register
    class LBBoundary(ScriptInterfaceHelper):
        """
        Creates a LB boundary.

        """

        _so_name = "LBBoundaries::LBBoundary"
        _so_bind_methods = ("get_force",)
