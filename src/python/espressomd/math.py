from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class CylTrafoParams(ScriptInterfaceHelper):
    _so_name = "CylTrafoParams"
