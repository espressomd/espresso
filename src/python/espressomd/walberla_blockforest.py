from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class WalberlaBlockForest(ScriptInterfaceHelper):
    """Interface to the Walberla BlockForest

    """
    _so_name = "walberla::WalberlaBlockForest"
    _so_creation_policy = "GLOBAL"

