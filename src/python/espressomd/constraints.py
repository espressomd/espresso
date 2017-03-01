from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper





class Constraints(ScriptInterfaceHelper):
    _so_name = "Constraints::Constraints"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Constraint):
                constraint = args[0]
            else:
                raise TypeError(
                    "Either a Constraint object or key-value pairs for the parameters of a Constraint object need to be passed.")
        else:
            constraint = Constraint(**kwargs)
        self.call_method("add", object=constraint)
        return constraint

    def remove(self, constraint):
        self.call_method("remove", constraint=constraint)

    def get(self):
        return self.call_method("get_elements")


class Constraint(ScriptInterfaceHelper):
    _so_name = "Constraints::Constraint"
