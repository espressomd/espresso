from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class Fene(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::Fene"
