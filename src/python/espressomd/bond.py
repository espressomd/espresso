from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class Fene(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::Fene"

@script_interface_register
class Harmonic(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::Harmonic"

@script_interface_register
class HarmonicDumbell(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::HarmonicDumbbell"

@script_interface_register
class BondedCoulomb(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::BondedCoulomb"
