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


@script_interface_register
class OverlapBondAngle(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::OverlapBondAngle"

@script_interface_register
class OverlapBondDihedral(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::OverlapBondDihedral"

@script_interface_register
class TabulatedBondAngle(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::TabulatedBondAngle"
    
@script_interface_register
class TabulatedBondDihedral(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::TabulatedBondDihedral"
    
@script_interface_register
class AngleHarmonic(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::AngleHarmonic"


@script_interface_register
class AngleCosine(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::AngleCosine"

    
@script_interface_register
class AngleCosSquare(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::AngleCosSquare"

@script_interface_register
class IbmTriel(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::IbmTriel"


@script_interface_register
class HydrogenBond(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::HydrogenBond"

@script_interface_register
class MembraneCollision(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::MembraneCollision"


@script_interface_register
class OifLocalForces(ScriptInterfaceHelper):
    """
    Docstring
    """

    _so_name = "Bond::OifLocalForces"
