from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper

class Constraints(ScriptInterfaceHelper):
    """
    List of active constraints. Add a :class:`espressomd.constraints.Constraint`
    to make it active in the system, or remove it to make it inactive.
    """

    _so_name = "Constraints::Constraints"

    def add(self, *args, **kwargs):
        """
        Add a constraint to the list.

        Parameters
        ----------
        Either an :class:`espressomd.constraints.Constraint`, or
        the parameters to construct one.
        """

        if len(args) == 1:
            if isinstance(args[0], ShapeBasedConstraint):
                print("option 1")
                constraint = args[0]
            else:
                raise TypeError(
                    "Either a Constraint object or key-value pairs for the parameters of a Constraint object need to be passed.")
        else:
            print("option 2")
            constraint = ShapeBasedConstraint(**kwargs)
        self.call_method("add", object=constraint)
        return constraint

    def remove(self, constraint):
        """
        Remove a constraint from the list.

        Parameters
        ----------
        constraint : :class:`espressomd.constraints.Constraint`
        """

        self.call_method("remove", constraint=constraint)


class Constraint(ScriptInterfaceHelper):
    """
    Attributes
    ----------
    only_positive : bool
      Act only in the direction of positive normal,
      only useful if penetrable is True.
    particle_type : int
      Interaction type of the constraint.
    penetrable : bool
      Whether particles are allowed to penetrate the
      constraint.
    shape : object
      One of the shapes from :mod:`espressomd.shapes`
    """
    _so_name = "Constraints::Constraint"

class ShapeBasedConstraint(ScriptInterfaceHelper):

    _so_name = "Constraints::ShapeBasedConstraint"
