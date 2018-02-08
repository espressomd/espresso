from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class VirtualSitesOff(ScriptInterfaceHelper):
    """Virtual sites implementation which does nothing (default)"""
    _so_name = "VirtualSites::VirtualSitesOff"



@script_interface_register
class VirtualSitesRelative(ScriptInterfaceHelper):
    """Virtual sites implementation placing virtual sites relative to other particles.
       See :ref:`Rigid arrangements of particles` for details. 

       Attributes
       ----------
       have_velocity : :obj:`bool`
           Determines whether the velocity of the virtual sites is calculated.
           This carries a performance cost.
       
       Attributes can be set on the instance or passed to the constructor as
       keyword arguments.

    """

       
    _so_name = "VirtualSites::VirtualSitesRelative"

@script_interface_register
class ActiveVirtualSitesHandle(ScriptInterfaceHelper):
    """Handle for the virtual sites implementation active in the core

       This should not be used directly.

       Attributes
       ----------
       implementation : instance of a virtual sites implementation

    """
    _so_name = "VirtualSites::ActiveVirtualSitesHandle"
