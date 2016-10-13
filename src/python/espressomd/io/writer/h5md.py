from ...script_interface import ScriptInterfaceHelper
import sys

class H5md(ScriptInterfaceHelper):
    _so_name = "ScriptInterface::Writer::H5mdScript"
    _so_bind_methods = ["init_file", "write", "close"]
