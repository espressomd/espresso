from __future__ import print_function, absolute_import
from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register


@script_interface_register
class EKContainer(ScriptObjectRegistry):

    """
    List of active constraints. Add a :class:`espressomd.constraints.Constraint`
    to make it active in the system, or remove it to make it inactive.

    """

    _so_name = "walberla::EKContainer"

    def add(self, *args):

        if len(args) == 1:
            if isinstance(args[0], EKSpecies):
                self.call_method("add", object=args[0])
            else:
                raise TypeError(
                    "EKSpecies object needs to be passed.")
        return args[0]

    def remove(self, ekspecies):
        self.call_method("remove", object=ekspecies)

    def clear(self):
        self.call_method("clear")


@script_interface_register
class EKSpecies(ScriptInterfaceHelper):
    """Interface to the Walberla EKSpecies
    """
    _so_name = "walberla::EKSpecies"
    _so_creation_policy = "GLOBAL"
