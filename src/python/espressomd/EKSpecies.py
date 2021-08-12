from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class EKSpecies(ScriptInterfaceHelper):
    """Interface to the Walberla EKSpecies
    """
    _so_name = "walberla::EKSpecies"
    _so_creation_policy = "GLOBAL"