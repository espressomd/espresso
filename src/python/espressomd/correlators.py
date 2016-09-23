
from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper

class Correlator(ScriptInterfaceHelper):
    _so_name = "Correlators::Correlator"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Correlator):
                Correlator = args[0]
            else:
                raise TypeError(
                    "Either a Correlator object or key-value pairs for the parameters of a Correlator object need to be passed.")
        else:
            Correlator = Correlator(**kwargs)
        self.call_method("add", object=Correlator)
        return Correlator

    def remove(self, Correlator):
        self.call_method("remove", object=Correlator)



class Correlator(ScriptInterfaceHelper):
    _so_name = "Correlators::Correlator"



