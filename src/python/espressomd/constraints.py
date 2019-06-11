# Copyright (C) 2010-2018 The ESPResSo project
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
      Interaction velocity of the boundary
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
    >>> # now create an un-penetrable shape-based constraint of type 0
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
    See also the examples on ForceField.

    Attributes
    ----------

    field_data: array_like :obj:`float`:
        The actual field please be aware that depending on the interpolation
        order additional points are used on the boundaries.

    grid_spacing: array_like :obj:`float`:
        The spacing of the grid points.

    """

    def __init__(self, field, **kwargs):
        shape, codim = self._unpack_dims(field)

        super(
            _Interpolated, self).__init__(_field_shape=shape, _field_codim=codim,
                                          _field_data=field.flatten(), **kwargs)

    @classmethod
    def required_dims(cls, box_size, grid_spacing):
        """Calculate the grid size and origin needed for specified box size and grid spacing.
           Returns the shape and origin (coordinates of [0][0][0]) needed.

        Arguments
        ---------
        box_size : array_like obj:`float`
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
        box_size : array_like obj:`float`
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
        """Returns an array of the coordinates of  the grid
        points required.

        Arguments
        ---------
        box_size : array_like obj:`float`
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
        their individual scaling factor before it is applied.

    """

    def __init__(self, field, **kwargs):
        super(ForceField, self).__init__(field, **kwargs)

    _codim = 3
    _so_name = "Constraints::ForceField"


@script_interface_register
class PotentialField(_Interpolated):

    """
    A generic tabulated force field that applies a per particle
    scaling factor. The forces are calculated numerically from
    the data by finite differences. The potential is interpolated
    from the provided data.

    Attributes
    ----------
    default_scale : :obj:`float`
        Scaling factor for particles that have no
        individual scaling factor.
    particle_scales: array_like (:obj:`int`, :obj:`float`)
        A list of tuples of ids and scaling factors. For
        particles in the list the interaction is scaled with
        their individual scaling factor before it is applied.

    """

    def __init__(self, field, **kwargs):
        super(PotentialField, self).__init__(field, **kwargs)

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

    resulting in the electric field E
    everywhere. (E.g. in a plate capacitor).
    The resulting force on the particles are then

      F = q * E

    where q is the charge of the particle.

    Attributes
    ----------
    E : array of :obj:`float`
        The electric field.

    phi0 : :obj:`float`
           The potential at the origin

    """

    def __init__(self, E, phi0=0):
        super(LinearElectricPotential, self).__init__(A=-E, b=phi0)

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

      E = E0 * sin(k * x + omega * t + phi)

    The resulting force on the particles are then

      F = q * E

    where q is the charge of the particle.
    This can be used to generate a homogeneous AC
    field by setting k to zero.

    Attributes
    ----------
    E0 : array of :obj:`float`
        The amplitude of the electric field.
    k  : array of :obj:`float`
        Wave vector of the wave
    omega : :obj:`float`
        Frequency of the wave
    phi : :obj:`float`
           Optional phase shift, defaults to 0.

    """

    _so_name = "Constraints::ElectricPlaneWave"

    def __init__(self, E0, k, omega, phi=0):
        super(ElectricPlaneWave, self).__init__(amplitude=E0,
                                                wave_vector=k,
                                                frequency=omega,
                                                phase=phi)

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

      F = -gamma * (u(r) - v)

    where v is the velocity of the particle.

    """

    def __init__(self, field, **kwargs):
        super(FlowField, self).__init__(field, **kwargs)

    _codim = 3
    _so_name = "Constraints::FlowField"


@script_interface_register
class HomogeneousFlowField(Constraint):

    """
    Viscous coupling to a flow field that is
    constant in space with the force

      F = -gamma * (u - v)

    where v is the velocity of the particle.

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
    Electric potential interpolated from
    provided data. The electric field E is
    calculated numerically from the potential,
    and the resulting force on the particles are

      F = q * E

    where q is the charge of the particle.


    """

    def __init__(self, field, **kwargs):
        super(ElectricPotential, self).__init__(field, **kwargs)

    _codim = 1
    _so_name = "Constraints::ElectricPotential"
