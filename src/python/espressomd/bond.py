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

@script_interface_register
class Quartic(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::Quartic"


@script_interface_register
class SubtLj(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::SubtLj"

@script_interface_register
class Umbrella(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::Umbrella"

    
@script_interface_register
class TabulatedBondLength(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::TabulatedBondLength"


@script_interface_register
class OverlapBondLength(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::OverlapBondLength"

@script_interface_register
class ThermalizedBond(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::ThermalizedBond"
