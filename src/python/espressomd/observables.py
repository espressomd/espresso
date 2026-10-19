#
# Copyright (C) 2010-2026 The ESPResSo project
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
#

import itertools
import numpy as np
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .math import CylindricalTransformationParameters


def _particles_to_ids(particles):
    """
    Convert particle selection to a list of integer particle ids.

    Parameters
    ----------
    particles : array_like
        Particle selection. Elements can be particle handles (objects exposing
        an ``id`` attribute), a particle slice (iterable of particle handles),
        or integer particle ids.

    Returns
    -------
    list of :obj:`int`
        Particle ids in the same order as provided.
    """
    if particles is None:
        raise TypeError("'particles' must not be None")

    # Single numeric id
    if isinstance(particles, (int, np.integer)):
        return [int(particles)]

    # Objects with attribute "id" (ParticleHandle OR ParticleSlice-like)
    if hasattr(particles, "id"):
        pid = particles.id

        # If id is scalar -> single particle
        if np.isscalar(pid) or isinstance(pid, (int, np.integer)):
            return [int(pid)]

        # If id is array-like -> slice selection
        try:
            arr = np.asarray(pid)
            if arr.ndim == 0:
                return [int(arr)]
            return [int(x) for x in arr.tolist()]
        except Exception:
            # Fall back to iterable handling below
            pass

    # Iterable of ids / handles
    try:
        iterator = iter(particles)
    except TypeError as e:
        raise TypeError(
            "'particles' must be an int, a particle handle, a particle slice, or an iterable of those"
        ) from e

    ids = []
    for item in iterator:
        if isinstance(item, (int, np.integer)):
            ids.append(int(item))
        elif hasattr(item, "id"):
            ids.append(int(item.id))
        else:
            raise TypeError(
                "Invalid element in 'particles': expected int or a particle handle (object with attribute 'id')"
            )
    return ids


@script_interface_register
class Observable(ScriptInterfaceHelper):
    """
    Base class for all observables.

    Methods
    -------
    shape()
        Get the shape of the numpy array returned by the observable.
    """
    _so_name = "Observables::Observable"
    _so_bind_methods = ("shape",)
    _so_creation_policy = "GLOBAL"

    # If defined in a subclass, maps public parameter name -> backend keyword name.
    # Example: {"particles": "ids"} or {"particles": "ids", "target_particles": "target_ids"}.
    _particle_param_map = None
    _optional_particle_params = ()

    def __init__(self, **kwargs):
        # Observables without particle selection: reject particle keywords
        particle_map = getattr(type(self), "_particle_param_map", None)
        if particle_map is None:
            forbidden = {"particles", "target_particles",
                         "particles1", "particles2"}
            used = forbidden.intersection(kwargs.keys())
            if used:
                raise TypeError(
                    f"{type(self).__name__} does not accept {sorted(used)}")
            super().__init__(**kwargs)
            return

        defaults = dict(
            getattr(type(self), "_particle_param_defaults", {}) or {})

        old_to_new = {
            "ids": "particles",
            "target_ids": "target_particles",
            "ids1": "particles1",
            "ids2": "particles2",
        }
        for old, new in old_to_new.items():
            if old in kwargs:
                raise TypeError(
                    f"Parameter '{old}' has been renamed to '{new}'")

        # Fill defaults for missing optional params (e.g., RDF particles2 -> [])
        for public_name, default_val in defaults.items():
            if public_name not in kwargs:
                kwargs[public_name] = default_val

        # Conversion for declared params
        for public_name, backend_name in type(self)._particle_param_map.items():
            if public_name in kwargs:
                kwargs[backend_name] = _particles_to_ids(
                    kwargs.pop(public_name))

        super().__init__(**kwargs)

    def calculate(self):
        return np.array(self.call_method("calculate")).reshape(self.shape())


class ProfileObservable(Observable):
    """
    Base class for histogram-based observables.
    """

    def bin_edges(self):
        """
        Returns
        -------
        :obj:`ndarray` of :obj:`float`
            Positions between the bins. If the histogram has dimensions
            ``(M,N,O)``, the bin edges have dimensions ``(M+1,N+1,O+1,3)``.
        """
        edges = self.call_method("edges")
        shape = list(map(len, edges)) + [len(edges)]
        return np.array(list(itertools.product(*edges))).reshape(shape)

    def bin_centers(self):
        """
        Returns
        -------
        :obj:`ndarray` of :obj:`float`
            Positions of the bins centers. If the histogram has dimensions
            ``(M,N,O)``, the bin centers have dimensions ``(M,N,O,3)``.
        """
        edges = self.call_method("edges")
        for i, edge in enumerate(edges):
            edges[i] = np.array(edge[:-1]) + (edge[1] - edge[0]) / 2
        shape = list(map(len, edges)) + [len(edges)]
        return np.array(list(itertools.product(*edges))).reshape(shape)


class CylindricalProfileObservable(ProfileObservable):
    """
    Base class for observables that work with cylinder coordinates
    """

    def __init__(self, transform_params=CylindricalTransformationParameters(),
                 **kwargs):
        # Provide default transformation parameters if not user-provided
        kwargs['transform_params'] = transform_params
        super().__init__(**kwargs)


@script_interface_register
class ComPosition(Observable):

    """Calculates the center of mass for particles with given ids.

    Note that virtual sites are not included since they do not have a meaningful mass.

    Output format: :math:`\\frac{1}{\\sum_i m_i} \\left( \\sum_i m_i r^x_i, \\sum_i m_i r^y_i, \\sum_i m_i r^z_i\\right)`

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ComPosition"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ComVelocity(Observable):

    """Calculates the center of mass velocity for particles with given ids.

    Note that virtual sites are not included since they do not have a meaningful mass.

    Output format: :math:`\\frac{1}{\\sum_i m_i} \\left( \\sum_i m_i v^x_i, \\sum_i m_i v^y_i, \\sum_i m_i v^z_i\\right)`

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ComVelocity"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class DensityProfile(ProfileObservable):

    """Calculates the particle density profile for particles with given ids.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider (inclusive).
    min_y : :obj:`float`
        Minimum ``y`` to consider (inclusive).
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_x : :obj:`float`
        Maximum ``x`` to consider (exclusive).
    max_y : :obj:`float`
        Maximum ``y`` to consider (exclusive).
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_x_bins``, ``n_y_bins``, ``n_z_bins``) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DensityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class DipoleMoment(Observable):

    """Calculates the electric dipole moment for particles with given ids.

    Output format: :math:`\\left(\\sum_i q_i r^x_i, \\sum_i q_i r^y_i, \\sum_i q_i r^z_i\\right)`

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DipoleMoment"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class FluxDensityProfile(ProfileObservable):

    """Calculates the particle flux density for particles with given ids.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider (inclusive).
    min_y : :obj:`float`
        Minimum ``y`` to consider (inclusive).
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_x : :obj:`float`
        Maximum ``x`` to consider (exclusive).
    max_y : :obj:`float`
        Maximum ``y`` to consider (exclusive).
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the x,
            y and z components of the flux density, respectively.

    """
    _so_name = "Observables::FluxDensityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ForceDensityProfile(ProfileObservable):

    """Calculates the force density profile for particles with given ids.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider (inclusive).
    min_y : :obj:`float`
        Minimum ``y`` to consider (inclusive).
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_x : :obj:`float`
        Maximum ``x`` to consider (exclusive).
    max_y : :obj:`float`
        Maximum ``y`` to consider (exclusive).
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the x,
            y and z components of the force, respectively.

    """
    _so_name = "Observables::ForceDensityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class LBVelocityProfile(ProfileObservable):

    """Calculates the LB fluid velocity profile.

    This observable samples the fluid in on a regular grid defined by the variables
    ``sampling_*``. Note that a small delta leads to a large number of sample
    points and carries a performance cost.

    Parameters
    ----------
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider (inclusive).
    min_y : :obj:`float`
        Minimum ``y`` to consider (inclusive).
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_x : :obj:`float`
        Maximum ``x`` to consider (exclusive).
    max_y : :obj:`float`
        Maximum ``y`` to consider (exclusive).
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).
    sampling_delta_x : :obj:`float`, default=1.0
        Spacing for the sampling grid in ``x``-direction.
    sampling_delta_y : :obj:`float`, default=1.0
        Spacing for the sampling grid in ``y``-direction.
    sampling_delta_z : :obj:`float`, default=1.0
        Spacing for the sampling grid in ``z``-direction.
    sampling_offset_x : :obj:`float`, default=0.0
        Offset for the sampling grid in ``x``-direction.
    sampling_offset_y : :obj:`float`, default=0.0
        Offset for the sampling grid in ``y``-direction.
    sampling_offset_z : :obj:`float`, default=0.0
        Offset for the sampling grid in ``z``-direction.
    allow_empty_bins : :obj:`bool`, default=False
        Whether or not to allow bins that will not be sampled at all.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the x,
            y and z components of the LB velocity, respectively.

    """
    _so_name = "Observables::LBVelocityProfile"


@script_interface_register
class LBFluidPressureTensor(Observable):

    """Calculates the average pressure tensor of the LB fluid for all nodes.

    Parameters
    ----------
    None

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::LBFluidPressureTensor"


@script_interface_register
class MagneticDipoleMoment(Observable):

    """Calculates the magnetic dipole moment for particles with given ids.

    Output format: :math:`\\left(\\sum_i \\mu^x_i, \\sum_i \\mu^y_i, \\sum_i \\mu^z_i\\right)`

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::MagneticDipoleMoment"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleAngularVelocities(Observable):

    """Calculates the angular velocity (omega) in the spaced-fixed frame of reference

    Output format: :math:`(\\omega^x_1,\\ \\omega^y_1,\\ \\omega^z_1),\\ (\\omega^x_2,\\ \\omega^y_2,\\ \\omega^z_2), \\dots,\\ (\\omega^x_n,\\ \\omega^y_n,\\ \\omega^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleAngularVelocities"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleBodyAngularVelocities(Observable):

    """Calculates the angular velocity (omega) in the particles'  body-fixed frame of reference.

    For each particle, the body-fixed frame of reference is obtained from the particle's
    orientation stored in the quaternions.

    Output format: :math:`(\\omega^x_1,\\ \\omega^y_1,\\ \\omega^z_1),\\ (\\omega^x_2,\\ \\omega^y_2,\\ \\omega^z_2), \\dots,\\ (\\omega^x_n,\\ \\omega^y_n,\\ \\omega^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleBodyAngularVelocities"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleBodyVelocities(Observable):

    """Calculates the particle velocity in the particles'  body-fixed frame of reference.

    For each particle, the body-fixed frame of reference is obtained from the particle's
    orientation stored in the quaternions.

    Output format: :math:`(v^x_1,\\ v^y_1,\\ v^z_1),\\ (v^x_2,\\ v^y_2,\\ v^z_2),\\ \\dots,\\ (v^x_n,\\ v^y_n,\\ v^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleBodyVelocities"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleForces(Observable):

    """Calculates the particle forces for particles with given ids.

    Output format: :math:`(f^x_1,\\ f^y_1,\\ f^z_1),\\ (f^x_2,\\ f^y_2,\\ f^z_2),\\ \\dots,\\ (f^x_n,\\ f^y_n,\\ f^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleForces"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticlePositions(Observable):

    """Calculates the particle positions for particles with given ids.

    Output format: :math:`(x_1,\\ y_1,\\ z_1),\\ (x_2,\\ y_2,\\ z_2),\\ \\dots,\\ (x_n,\\ y_n,\\ z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticlePositions"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleVelocities(Observable):

    """Calculates the particle velocities for particles with given ids.

    Output format: :math:`(v^x_1,\\ v^y_1,\\ v^z_1),\\ (v^x_2,\\ v^y_2,\\ v^z_2),\\ \\dots,\\ (v^x_n,\\ v^y_n,\\ v^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleVelocities"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleDirectors(Observable):

    """Calculates the particle directors for particles with given ids.

    Output format: :math:`(d^x_1,\\ d^y_1,\\ d^z_1),\\ (d^x_2,\\ d^y_2,\\ d^z_2),\\ \\dots,\\ (d^x_n,\\ d^y_n,\\ d^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleDirectors"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleDipoleFields(Observable):

    """Calculates the particle dipole fields for particles with given ids.

    Output format: :math:`(h^x_1,\\ h^y_1,\\ h^z_1),\\ (h^x_2,\\ h^y_2,\\ h^z_2),\\ \\dots,\\ (h^x_n,\\ h^y_n,\\ h^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleDipoleFields"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class ParticleDistances(Observable):

    """Calculates the distances between particles with given ids along a
    polymer chain.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N - 1,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleDistances"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class PairwiseDistances(Observable):
    """
    Calculates the distance matrix between two sets of particles.
    Duplicate entries in each set are only counted once.
    The calculation yields a flattened triangular matrix.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The first set of ids of particles.

    target_particles : array_like of :obj:`int`
        The second set of (target) ids of particles.
        In case of overlap with the first set,
        self-interactions are removed from the result.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N * M / 2,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::PairwiseDistances"
    _particle_param_map = {"particles": "ids",
                           "target_particles": "target_ids"}


@script_interface_register
class TotalForce(Observable):

    """Calculates the total force on particles with given ids.

    Note that virtual sites are not included since forces on them do not enter the equation of motion directly.

    Output format: :math:`\\left(\\sum_i f^x_i, \\sum_i f^y_i, \\sum_i f^z_i\\right)`

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::TotalForce"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class BondAngles(Observable):

    """Calculates the angles between bonds of particles with given ids along a
    polymer chain.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N - 2,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::BondAngles"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CosPersistenceAngles(Observable):

    """Calculates the cosine of mutual bond angles for chained particles with given ids.

    The *i*-th  value of the result contains the cosine of the angle between bonds that
    are separated by *i* bonds. The values are averaged over the chain.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N - 2,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::CosPersistenceAngles"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class BondDihedrals(Observable):

    """Calculates the dihedrals between particles with given ids along a
    polymer chain.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (N - 3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::BondDihedrals"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class Energy(Observable):

    """Calculates the total energy.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        :obj:`float`

    """
    _so_name = "Observables::Energy"


@script_interface_register
class Pressure(Observable):

    """Calculates the total scalar pressure.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        :obj:`float`

    """
    _so_name = "Observables::Pressure"


@script_interface_register
class PressureTensor(Observable):

    """Calculates the total pressure tensor.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::PressureTensor"


@script_interface_register
class DPDPressure(Observable):

    """Calculates the non-equilibrium contribution of the DPD interaction
    to the pressure tensor.

    Parameters
    ----------
    None

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DPDPressure"


@script_interface_register
class CylindricalDensityProfile(CylindricalProfileObservable):

    """Calculates the particle density in cylindrical coordinates.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::CylindricalDensityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CylindricalFluxDensityProfile(CylindricalProfileObservable):

    """Calculates the particle flux density in cylindrical coordinates.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the
            radial distance, azimuth and axial coordinate of the particle
            flux density field, respectively.

    """
    _so_name = "Observables::CylindricalFluxDensityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CylindricalLBFluxDensityProfileAtParticlePositions(
        CylindricalProfileObservable):

    """Calculates the LB fluid flux density at the particle positions in
    cylindrical coordinates.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the
            radial distance, azimuth and axial coordinate of the LB flux
            density field, respectively.

    """
    _so_name = "Observables::CylindricalLBFluxDensityProfileAtParticlePositions"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CylindricalLBVelocityProfileAtParticlePositions(
        CylindricalProfileObservable):

    """Calculates the LB fluid velocity at the particle positions in
    cylindrical coordinates.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the
            radial distance, azimuth and axial coordinate of the LB velocity
            field, respectively.

    """
    _so_name = "Observables::CylindricalLBVelocityProfileAtParticlePositions"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CylindricalVelocityProfile(CylindricalProfileObservable):

    """Calculates the particle velocity profile in cylindrical coordinates.

    Parameters
    ----------
    particles : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the
            radial distance, azimuth and axial coordinate of the particle
            velocity field, respectively.

    """
    _so_name = "Observables::CylindricalVelocityProfile"
    _particle_param_map = {"particles": "ids"}


@script_interface_register
class CylindricalLBVelocityProfile(CylindricalProfileObservable):

    """Calculates the LB fluid velocity profile in cylindrical coordinates.

    This observable samples the fluid in on a regular grid defined by variable
    ``sampling_density``. Note that a small delta leads to a large number of
    sample points and carries a performance cost.

    Parameters
    ----------
    transform_params : :class:`espressomd.math.CylindricalTransformationParameters`, optional
        Parameters of the cylinder transformation. Defaults to the default of :class:`espressomd.math.CylindricalTransformationParameters`
    n_r_bins : :obj:`int`, default = 1
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`, default = 1
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`, default = 1
        Number of bins in ``z`` direction.
    min_r : :obj:`float`, default = 0
        Minimum ``r`` to consider (inclusive).
    min_phi : :obj:`float`, default = :math:`-\\pi`
        Minimum ``phi`` to consider (inclusive). Must be in :math:`[-\\pi,\\pi)`.
    min_z : :obj:`float`
        Minimum ``z`` to consider (inclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).
    max_phi : :obj:`float`, default = :math:`\\pi`
        Maximum ``phi`` to consider (exclusive). Must be in :math:`(-\\pi,\\pi]`.
    max_z : :obj:`float`
        Maximum ``z`` to consider (exclusive).
    sampling_density : :obj:`float`
        Samples per unit volume for the LB velocity interpolation.

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
            The fourth dimension of the array stores the histogram for the
            radial distance, azimuth and axial coordinate of the LB velocity
            field, respectively.

    """
    _so_name = "Observables::CylindricalLBVelocityProfile"


@script_interface_register
class RDF(Observable):

    """Calculates a radial distribution function.
    The result is normalized by the bulk concentration.

    Parameters
    ----------
    particles1 : array_like of :obj:`int`
        The ids of (existing) particles to calculate the distance from.
    particles2 : array_like of :obj:`int`, optional
        The ids of (existing) particles to calculate the distance to.
        If not provided, use ``particles1``.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider (exclusive).
    max_r : :obj:`float`
        Maximum ``r`` to consider (exclusive).

    Methods
    -------
    calculate()
        Run the observable.

        Returns
        -------
        (``n_r_bins``,) :obj:`ndarray` of :obj:`float`
            The RDF.

    """
    _so_name = "Observables::RDF"
    _particle_param_map = {"particles1": "ids1", "particles2": "ids2"}

    def __init__(self, **kwargs):
        # Preserve prior behavior: if second set not provided, backend uses ids1.
        if "particles2" not in kwargs:
            kwargs["particles2"] = []
        super().__init__(**kwargs)

    def bin_centers(self):
        bin_width = (self.max_r - self.min_r) / self.n_r_bins
        return self.min_r + (np.arange(self.n_r_bins) + 0.5) * bin_width
