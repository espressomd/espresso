from .script_interface import ScriptInterfaceHelper, script_interface_register

@script_interface_register
class VirtualSitesOff(ScriptInterfaceHelper):
    _so_name = "VirtualSites::VirtualSitesOff"

@script_interface_register
class VirtualSitesRelative(ScriptInterfaceHelper):
    _so_name = "VirtualSites::VirtualSitesRelative"

@script_interface_register
class ActiveVirtualSitesHandle(ScriptInterfaceHelper):
    _so_name = "VirtualSites::ActiveVirtualSitesHandle"
