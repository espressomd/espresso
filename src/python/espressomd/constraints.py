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
        the parameters to construct an :class:`espressomd.constraints.ShapeBasedConstraint`.
        """

        if len(args) == 1:
            if isinstance(args[0], Constraint):
                constraint = args[0]
            else:
                raise TypeError(
                    "Either a Constraint object or key-value pairs for the parameters of a ShapeBasedConstraint object need to be passed.")
        else:
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
    Base class for constraints. A constraint provides a force and
    an energy contribution for a single particle.
    """

    _so_name = "Constraints::Constraint"

class ShapeBasedConstraint(Constraint):
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

    See Also
    shapes : shape objects that define mathematical surfaces
    ----------

    Example
    ----------
    >>> import espressomd
    >>> from espressomd import shapes
    >>> system = espressomd.System()
    >>>
    >>> # create first a shape-object to define the constraint surface
    >>> spherical_cavity = shapes.Sphere(center=[5,5,5], radius=5.0, direction=-1.0)
    >>>
    >>> # now create an un-penetrable shape-based contraint of type 0
    >>> spherical_constraint = system.constraints.add(particle_type=0, penetrable=0, shape=spherical_cavity)
    >>>
    >>> #place a trapped particle inside this sphere
    >>> system.part.add(id=0, pos=[5,5,5], type=1)
    >>>

    
    """

    _so_name = "Constraints::ShapeBasedConstraint"
