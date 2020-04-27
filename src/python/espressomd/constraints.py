# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
import numpy as np
from itertools import product


@script_interface_register
class Constraints(ScriptObjectRegistry):

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
        constraint: :class:`espressomd.constraints.Constraint`
            Either a constraint object...
        \*\*kwargs : any
            ... or parameters to construct an
            :class:`espressomd.constraints.ShapeBasedConstraint`

        Returns
        ----------
        constraint : :class:`espressomd.constraints.Constraint`
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
        constraint : :obj:`espressomd.constraints.Constraint`

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
    only_positive : :obj:`bool`
        Act only in the direction of positive normal,
        only useful if penetrable is ``True``.
    particle_type : :obj:`int`
        Interaction type of the constraint.
    particle_velocity : array_like of :obj:`float`
        Interaction velocity of the boundary
    penetrable : :obj:`bool`
        Whether particles are allowed to penetrate the constraint.
    shape : :class:`espressomd.shapes.Shape`
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
    >>> # now create an un-penetrable shape-based constraint of type 0
    >>> spherical_constraint = system.constraints.add(particle_type=0, penetrable=False, shape=spherical_cavity)
    >>>
    >>> # place a trapped particle inside this sphere
    >>> system.part.add(id=0, pos=[5, 5, 5], type=1)

    """

    _so_name = "Constraints::ShapeBasedConstraint"

    def min_dist(self):
        """
        Calculates the minimum distance to all interacting particles.

        Returns
        ----------
        :obj:`float` :
            The minimum distance
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
        ...     epsilon=1, sigma=1,
        ...     cutoff=2**(1. / 6), shift="auto")
        >>>
        >>> floor = system.constraints.add(shape=shapes.Wall(normal=[0, 0, 1], dist=0.0),
        ...    particle_type=0, penetrable=False, only_positive=False)
        >>>
        >>> system.part.add(id=0, pos=[0,0,1.5], type=0, ext_force=[0, 0, -.1])
        >>> # print the particle position as it falls
        >>> # and print the force it applies on the floor
        >>> for t in range(10):
        ...     system.integrator.run(100)
        ...     print(system.part[0].pos, floor.total_force())

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
    H : (3,) array_like of :obj:`float`
        Magnetic field vector. Describes both field direction and
        strength of the magnetic field (via length of the vector).

    """

    _so_name = "Constraints::HomogeneousMagneticField"


class _Interpolated(Constraint):

    """
    Tabulated field data.
    The actual field value is calculated by linear
    interpolation (force fields) or gradient linear
    interpolation.

    The data has to have one point of halo in each direction,
    and is shifted by half a grid spacing in the +xyz direction,
    so that the element (0,0,0) has coordinates -0.5 * grid_spacing.
    The number of points has to be such that the data spans the whole
    box, e.g. the most up right back point has to be at least at
    box + 0.5 * grid_spacing. There are convenience functions on this
    class that can calculate the required grid dimensions and the coordinates.

    Arguments
    ----------
    field : (M, N, O, P) array_like of :obj:`float`
        The actual field on a grid of size (M, N, O) with dimension P.
    grid_spacing : (3,) array_like of :obj:`float`
        Spacing of the grid points.

    Attributes
    ----------

    field : (M, N, O, P) array_like of :obj:`float`
        The actual field on a grid of size (M, N, O) with dimension P.
        Please be aware that depending on the interpolation
        order additional points are used on the boundaries.

    grid_spacing : array_like of :obj:`float`
        Spacing of the grid points.

    origin : (3,) array_like of :obj:`float`
        Coordinates of the grid origin.

    """

    def __init__(self, **kwargs):
        if "oid" not in kwargs:
            field = kwargs.pop("field")
            shape, codim = self._unpack_dims(field)
            super().__init__(_field_shape=shape, _field_codim=codim,
                             _field_data=field.flatten(), **kwargs)
        else:
            super().__init__(**kwargs)

    @classmethod
    def required_dims(cls, box_size, grid_spacing):
        """
        Calculate the grid size and origin needed for specified box size and
        grid spacing. Returns the shape and origin (coordinates of [0][0][0])
        needed.

        Arguments
        ---------
        box_size : (3,) array_like of obj:`float`
            The box the field should be used.

        grid_spacing : array_like obj:`float`
            The desired grid spacing.

        """

        shape = np.array(np.ceil(box_size / grid_spacing), dtype=int) + 2
        origin = -0.5 * grid_spacing
        return shape, origin

    @classmethod
    def field_from_fn(cls, box_size, grid_spacing, f, codim=None):
        """Generate field data for a desired box size and grid spacing
        by evaluating a function at the coordinates.

        Arguments
        ---------
        box_size : (3,) array_like of obj:`float`
            The box the field should be used.

        grid_spacing : array_like obj:`float`
            The desired grid spacing.

        f : callable
           A function that is called with the coordinates of every
           grid point to populate the grid.

        """

        shape, origin = cls.required_dims(box_size, grid_spacing)

        if not codim:
            codim = cls._codim

        field = np.zeros((shape[0], shape[1], shape[2], codim))

        for i in product(*map(range, shape)):
            x = origin + np.array(i) * grid_spacing
            field[i] = f(x)

        return field

    @classmethod
    def field_coordinates(cls, box_size, grid_spacing):
        """Returns an array of the coordinates of the grid points required.

        Arguments
        ---------
        box_size : (3,) array_like of obj:`float`
            The box the field should be used.

        grid_spacing : array_like obj:`float`
            The desired grid spacing.
        """

        return cls.field_from_fn(box_size, grid_spacing, lambda x: x, 3)

    def _unpack_dims(self, a):
        s = a.shape
        shape = s[:3]
        codim = s[3]

        return (shape, codim)

    @property
    def field(self):
        shape = self._field_shape
        return np.reshape(self._field_data,
                          (shape[0], shape[1], shape[2], self._field_codim))


@script_interface_register
class ForceField(_Interpolated):

    """
    A generic tabulated force field that applies a per-particle scaling factor.

    Arguments
    ----------
    field : (M, N, O, 3) array_like of :obj:`float`
        Forcefield amplitude on a grid of size (M, N, O).
    grid_spacing : (3,) array_like of :obj:`float`
        Spacing of the grid points.
    default_scale : :obj:`float`
        Scaling factor for particles that have no individual scaling factor.
    particle_scales : array_like of (:obj:`int`, :obj:`float`)
        A list of tuples of ids and scaling factors. For
        particles in the list the interaction is scaled with
        their individual scaling factor before it is applied.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    _codim = 3
    _so_name = "Constraints::ForceField"


@script_interface_register
class PotentialField(_Interpolated):

    """
    A generic tabulated force field that applies a per-particle
    scaling factor. The forces are calculated numerically from
    the data by finite differences. The potential is interpolated
    from the provided data.

    Arguments
    ----------
    field : (M, N, O, 1) array_like of :obj:`float`
        Potential on a grid of size (M, N, O).
    grid_spacing : (3,) array_like of :obj:`float`
        Spacing of the grid points.
    default_scale : :obj:`float`
        Scaling factor for particles that have no individual scaling factor.
    particle_scales : array_like (:obj:`int`, :obj:`float`)
        A list of tuples of ids and scaling factors. For
        particles in the list the interaction is scaled with
        their individual scaling factor before it is applied.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    _codim = 1
    _so_name = "Constraints::PotentialField"


@script_interface_register
class Gravity(Constraint):

    """
    Gravity force

    :math:`F = m \\cdot g`

    Arguments
    ----------
    g : (3,) array_like of :obj:`float`
        The gravitational acceleration.

    """

    def __init__(self, **kwargs):
        if "oid" not in kwargs:
            kwargs["value"] = kwargs.pop("g")
        super().__init__(**kwargs)

    @property
    def g(self):
        return self.value

    _so_name = "Constraints::Gravity"


@script_interface_register
class LinearElectricPotential(Constraint):

    """
    Electric potential of the form

    :math:`\\phi = -E \\cdot x + \\phi_0`,

    resulting in the electric field E
    everywhere. (E.g. in a plate capacitor).
    The resulting force on the particles are then

    :math:`F = q \\cdot E`

    where :math:`q` is the charge of the particle.

    Arguments
    ----------
    E : array_like of :obj:`float`
        The electric field.

    phi0 : :obj:`float`
        The potential at the origin

    """

    def __init__(self, phi0=0, **kwargs):
        if "oid" not in kwargs:
            kwargs["A"] = -np.array(kwargs.pop("E"))
            kwargs["b"] = phi0
        super().__init__(**kwargs)

    @property
    def E(self):
        return -np.array(self.A)

    @property
    def phi0(self):
        return np.array(self.b)

    _so_name = "Constraints::LinearElectricPotential"


@script_interface_register
class ElectricPlaneWave(Constraint):

    """
    Electric field of the form

    :math:`E = E0 \\cdot \\sin(k \\cdot x + \\omega \\cdot t + \\phi)`

    The resulting force on the particles are then

    :math:`F = q \\cdot E`

    where :math:`q` is the charge of the particle.
    This can be used to generate a homogeneous AC
    field by setting k to zero.

    Arguments
    ----------
    E0 : array_like of :obj:`float`
        Amplitude of the electric field.
    k  : array_like of :obj:`float`
        Wave vector of the wave
    omega : :obj:`float`
        Frequency of the wave
    phi : :obj:`float`, optional
        Phase shift

    """

    _so_name = "Constraints::ElectricPlaneWave"

    def __init__(self, phi=0, **kwargs):
        if "oid" not in kwargs:
            kwargs["amplitude"] = kwargs.pop("E0")
            kwargs["wave_vector"] = kwargs.pop("k")
            kwargs["frequency"] = kwargs.pop("omega")
            kwargs["phase"] = phi
        super().__init__(**kwargs)

    @property
    def E0(self):
        return np.array(self.amplitude)

    @property
    def k(self):
        return np.array(self.wave_vector)

    @property
    def omega(self):
        return self.frequency

    @property
    def phi(self):
        return self.phase


@script_interface_register
class FlowField(_Interpolated):

    """
    Viscous coupling to a flow field that is
    interpolated from tabulated data like

    :math:`F = -\\gamma \\cdot \\left( u(r) - v \\right)`

    where :math:`v` is the velocity of the particle.

    Arguments
    ----------
    field : (M, N, O, 3) array_like of :obj:`float`
        Field velocity on a grid of size (M, N, O)
    grid_spacing : (3,) array_like of :obj:`float`
        Spacing of the grid points.
    gamma : :obj:`float`
        Coupling constant

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    _codim = 3
    _so_name = "Constraints::FlowField"


@script_interface_register
class HomogeneousFlowField(Constraint):

    """
    Viscous coupling to a flow field that is
    constant in space with the force

    :math:`F = -\\gamma \\cdot (u - v)`

    where :math:`v` is the velocity of the particle.

    Attributes
    ----------
    gamma : :obj:`float`
        Coupling constant

    """

    def __init__(self, **kwargs):
        if "oid" not in kwargs:
            kwargs["value"] = kwargs.pop("u")
        super().__init__(**kwargs)

    @property
    def u(self):
        """
        Field velocity ((3,) array_like of :obj:`float`).
        """
        return self.value

    _so_name = "Constraints::HomogeneousFlowField"


@script_interface_register
class ElectricPotential(_Interpolated):

    """
    Electric potential interpolated from
    provided data. The electric field E is
    calculated numerically from the potential,
    and the resulting force on the particles are

    :math:`F = q \\cdot E`

    where :math:`q` is the charge of the particle.

    Arguments
    ----------
    field : (M, N, O, 1) array_like of :obj:`float`
        Potential on a grid of size (M, N, O)
    grid_spacing : (3,) array_like of :obj:`float`
        Spacing of the grid points.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    _codim = 1
    _so_name = "Constraints::ElectricPotential"
