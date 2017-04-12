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

        if isinstance(args[0], Constraint):
            constraint = args[0]
        else:
            constraint = Constraint(**kwargs)
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
    ext_electric_field : double
      External electric field magnitude along the
      constraint shape surface normal direction.
      This feature is currently precisely implemented
      for the :mod:`espressomd.shapes.Wall` only.
    ext_magn_field : double
      External magnetic field magnitude along the
      constraint shape surface normal direction.
      This feature is currently precisely implemented
      for the :mod:`espressomd.shapes.Wall` only.
    """
    _so_name = "Constraints::Constraint"
