from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
from espressomd.utils import is_valid_type
import numpy as np
from itertools import product

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

@script_interface_register
class HomogeneousMagneticField(Constraint):
    """
    Attributes
    ----------
    H : array of :obj:`float`
        Magnetic field vector. Describes both field direction and
        strength of the magnetic field (via length of the vector).

    """

    _so_name = "Constraints::HomogeneousMagneticField"

class _Interpolated(Constraint):
    """
    Tabulated field data.
    The actual field value is calculated by cardinal b-spline
    intepolation (force fields) or gradient cardinal b-splins
    interpolation.

    The provided field data grid has to be larger than the box
    by a margin depending on the interpolation order:
    In each dimension the provided grid has to be
    order//2 * grid_spacing larger than the box.
    For example if the box is 10 in x direction, and the interpolation
    order is 2 and a grid spacing of .1 is to be used, (2//2) = 1 extra
    point is needed on each side, and the grid spans the range
    [-1 * 0.1, 10 + 1 * 0.1]. Please be aware that the periodicity is not
    taken into account automatically, so this is also true for periodic
    directions.

    Attributes
    ----------
    interpolation_order: array_like :obj:`int`
        The order of the b-splines top be used. This is equivalent
        to the number of points to be used in each direction. E.g.
        order 2 corresponds to linear interpolation, 3 to quadratic
        and so on.

    field_data: array_like :obj:`float`:
        The actual field please be aware that depending on the interpolation
        order additional points are used on the boundaries.

    grid_spacing: array_like :obj:`float`:
        The spacing of the grid points.

    """

    def __init__(self, field, **kwargs):
        shape, codim = self._unpack_dims(field)

        super(_Interpolated, self).__init__(_field_shape=shape, _field_codim=codim,
                                         _field_data=field.flatten(), **kwargs)

    @classmethod
    def required_dims(cls, box_size, grid_spacing, order):
        """Calculate the grid size needed for specified box size, grid spacing and order.
        """
        halo_points = (order) // 2
        shape= np.array(np.ceil(box_size/grid_spacing), dtype=int) + 2 * halo_points
        origin = np.array(-(halo_points + 0.5)*grid_spacing)

        return shape, origin

    @classmethod
    def field_from_fn(cls, box_size, grid_spacing, order, f, codim=None):
        shape, origin = cls.required_dims(box_size, grid_spacing, order)

        if not codim:
            codim = cls._codim

        field = np.zeros((shape[0], shape[1], shape[2], codim))

        for i in product(*map(range,shape)):
            x = origin + np.array(i)*grid_spacing
            field[i] = f(x)

        return field

    @classmethod
    def field_coordinates(cls, box_size, grid_spacing, order):
        return cls.field_from_fn(box_size, grid_spacing, order, lambda x: x, 3)

    def _unpack_dims(self, a):
        s = a.shape
        shape = s[:3]
        codim = s[3]

        return (shape, codim)

    @property
    def field(self):
        shape = self._field_shape
        return np.reshape(self._field_data, (shape[0], shape[1], shape[2], self._field_codim))

@script_interface_register
class ForceField(_Interpolated):
    """
    A generic tabulated force field that applies a per particle
    scaling factor.

    Attributes
    ----------
    default_scale : :obj:`float`
        Scaling factor for particles that have no
        individual scaling factor.
    particle_scales: array_like (:obj:`int`, :obj:`float`)
        A list of tuples of ids and scaling factors. For
        particles in the list the interaction is scaled with
        their individual scaling factor befor it is applied.

    """

    def __init__(self, field, **kwargs):
        super(ForceField, self).__init__(field, **kwargs)

    _codim = 3
    _so_name = "Constraints::ForceField"


@script_interface_register
class PotentialField(_Interpolated):
    """
    A generic tabulated force field that applies a per particle
    scaling factor.

    Attributes
    ----------
    default_scale : :obj:`float`
        Scaling factor for particles that have no
        individual scaling factor.
    particle_scales: array_like (:obj:`int`, :obj:`float`)
        A list of tuples of ids and scaling factors. For
        particles in the list the interaction is scaled with
        their individual scaling factor befor it is applied.

    """

    def __init__(self, field, **kwargs):
        super(PotentialField, self).__init__(**kwargs)

    _codim = 1
    _so_name = "Constraints::PotentialField"


@script_interface_register
class Gravity(Constraint):
    """
    Gravity force
      F = m * g

    Attributes
    ----------
    g : array of :obj:`float`
        The gravitational acceleration.

    """

    def __init__(self, g):
        super(Gravity, self).__init__(value=g)

    @property
    def g(self):
        return self.value

    _so_name = "Constraints::Gravity"


@script_interface_register
class LinearElectricPotential(Constraint):
    """
    Electric potential of the form

      phi = -E * x + phi0,

    resulting in the electic field E
    everywhere. (E.g. in a plate capacitor).


    Attributes
    ----------
    E : array of :obj:`float`
        The electric field.

    phi0 : :obj:`float`
           The potential at the origin

    """

    def __init__(self, E, phi0 = 0):
        super(LinearElectricPotential, self).__init__(A=-E, b=phi0)

    @property
    def E(self):
        return -np.array(self.A)

    @property
    def phi0(self):
        return np.array(self.b)

    _so_name = "Constraints::LinearElectricPotential"


@script_interface_register
class FlowField(_Interpolated):
    """
    Viscous coupling to a flow field that is
    interpolated from tabulated data like

      F = -gamma * (u(r) - v)

    wher v is the velocity of the particle.

    """

    def __init__(self, field, **kwargs):
        super(FlowField, self).__init__(**kwargs)

    _codim = 3
    _so_name = "Constraints::FlowField"


@script_interface_register
class HomogeneousFlowField(Constraint):
    """
    Viscous coupling to a flow field that is
    constant in space with the force

      F = -gamma * (u - v)

    wher v is the velocity of the particle.

    Attributes
    ----------
    gamma : :obj:`float`
        The coupling constant
    u : array_like :obj:`float`
        The velocity of the field.

    """

    def __init__(self, u, gamma):
        super(HomogeneousFlowField, self).__init__(value=u, gamma=gamma)

    @property
    def u(self):
        return self.value

    _so_name = "Constraints::HomogeneousFlowField"

@script_interface_register
class ElectricPotential(_Interpolated):
    """
    Electric potential.

    """

    def __init__(self, field, **kwargs):
        super(ElectricPotential, self).__init__(**kwargs)

    _codim = 1
    _so_name = "Constraints::ElectricPotential"
