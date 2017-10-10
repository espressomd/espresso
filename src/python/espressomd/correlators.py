
from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper,script_interface_register
import numpy as np

@script_interface_register
class Correlator(ScriptInterfaceHelper):
    _so_name = "Correlators::Correlator"
    _so_bind_methods = ("update","auto_update","finalize", "dim_corr","n_results","hierarchy_depth")
    _so_creation_policy = "LOCAL"

    def result(self):
        res=np.array(self.call_method("get_correlation"))
        res=res.reshape((self.n_results(),2+self.dim_corr()))
        return res

@script_interface_register
class AutoUpdateCorrelators(ScriptInterfaceHelper):
    _so_name = "Correlators::AutoUpdateCorrelators"
    _so_creation_policy = "LOCAL"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Correlator):
                correlator = args[0]
            else:
                raise TypeError(
                    "Either a Correlator object or key-value pairs for the parameters of a Correlator object need to be passed.")
        else:
            correlator = Correlator(**kwargs)
        self.call_method("add", object=correlator)
        return correlator

    def remove(self, Correlator):
        self.call_method("remove", object=Correlator)
