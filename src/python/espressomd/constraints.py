from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
from espressomd.utils import is_valid_type


@script_interface_register
class Constraints(ScriptInterfaceHelper):
    """
    List of active constraints. Add a :class:`espressomd.constraints.Constraint`
    to make it active in the system, or remove it to make it inactive.

    """

    _so_name = "Constraints::Constraints"

    def __getitem__(self, key):
        return self.call_method("get_elements")[key]

    def __iter__(self):
        elements = self.call_method("get_elements")
        for e in elements:
            yield e

    def add(self, *args, **kwargs):
        """
        Add a constraint to the list.

        Parameters
        ----------
        Either an instance of :class:`espressomd.constraints.Constraint`, or
        the parameters to construct an :class:`espressomd.constraints.ShapeBasedConstraint`.

        Returns
        ----------
        constraint : Instance of :class:`espressomd.constraints.Constraint`
                     The added constraint

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
        constraint : Instance of :class:`espressomd.constraints.Constraint`

        """

        self.call_method("remove", object=constraint)

    def clear(self):
        """
        Remove all constraints.

        """
        self.call_method("clear")

class Constraint(ScriptInterfaceHelper):
    """
    Base class for constraints. A constraint provides a force and
    an energy contribution for a single particle.

    """

    _so_name = "Constraints::Constraint"


@script_interface_register
class ShapeBasedConstraint(Constraint):
    """

    Attributes
    ----------
    only_positive : bool
      Act only in the direction of positive normal,
      only useful if penetrable is True.
    particle_type : int
      Interaction type of the constraint.
    particle_velocity : array of :obj:`float`
      Interaction velocity of the boudary
    penetrable : bool
      Whether particles are allowed to penetrate the
      constraint.
    shape : object
      One of the shapes from :mod:`espressomd.shapes`

    See Also
    ----------
    espressomd.shapes : shape module that define mathematical surfaces

    Examples
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


    def min_dist(self):
        """
        Calculates the minimum distance to all interacting particles.

        Returns
        ----------
        :obj:float: The minimum distance
        """
        return self.call_method("min_dist", object=self)


    def total_force(self):
        """
        Get total force acting on this constraint.

        Examples
        ----------
        >>> import espressomd
        >>> from espressomd import shapes
        >>> system = espressomd.System()
        >>>
        >>> system.time_step = 0.01
        >>> system.box_l = [50, 50, 50]
        >>> system.thermostat.set_langevin(kT=0.0, gamma=1.0)
        >>> system.cell_system.set_n_square(use_verlet_lists=False)
        >>> system.non_bonded_inter[0, 0].lennard_jones.set_params(
        >>>     epsilon=1, sigma=1,
        >>>     cutoff=2**(1. / 6), shift="auto")
        >>>
        >>>
        >>> floor = system.constraints.add(shape=shapes.Wall(normal=[0, 0, 1], dist=0.0),
        >>>    particle_type=0, penetrable=0, only_positive=0)

        >>> system.part.add(id=0, pos=[0,0,1.5], type=0, ext_force=[0,0,-.1])
        >>> # print the particle position as it falls
        >>> # and print the force it applies on the floor
        >>> for t in range(10):
        >>>     system.integrator.run(100)
        >>>     print(system.part[0].pos, floor.total_force())

        """
        return self.call_method("total_force", constraint=self)
        
    def total_normal_force(self):
        """
        Get the total summed normal force acting on this constraint.
        
        """
        return self.call_method("total_normal_force", constraint=self)

class HomogeneousMagneticField(Constraint):
    """
    Attributes
    ----------
    H : array of :obj:`float`
        Magnetic field vector. Describes both field direction and
        strength of the magnetic field (via length of the vector).

    """

    _so_name = "Constraints::HomogeneousMagneticField"
