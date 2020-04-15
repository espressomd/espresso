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
import itertools
import numpy as np
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class Observable(ScriptInterfaceHelper):
    """
    Base class for all observables.

    Methods
    -------
    shape()
        Return the shape of the observable.
    """
    _so_name = "Observables::Observable"
    _so_bind_methods = ("shape",)
    _so_creation_policy = "LOCAL"

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


@script_interface_register
class ComPosition(Observable):

    """Calculates the center of mass for particles with given ids.

    Note that virtual sites are not included since they do not have a meaningful mass.

    Output format: :math:`\\frac{1}{\\sum_i m_i} \\left( \\sum_i m_i r^x_i, \\sum_i m_i r^y_i, \\sum_i m_i r^z_i\\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ComPosition"


@script_interface_register
class ComVelocity(Observable):

    """Calculates the center of mass velocity for particles with given ids.

    Note that virtual sites are not included since they do not have a meaningful mass.

    Output format: :math:`\\frac{1}{\\sum_i m_i} \\left( \\sum_i m_i v^x_i, \\sum_i m_i v^y_i, \\sum_i m_i v^z_i\\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ComVelocity"


@script_interface_register
class Current(Observable):

    """Calculates the electric current for particles with given ids.

    Output format: :math:`\\left(\\sum_i q_i v^x_i, \\sum_i q_i v^y_i, \\sum_i q_i v^z_i, \\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::Current"


@script_interface_register
class DensityProfile(ProfileObservable):

    """Calculates the particle density profile for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider.
    min_y : :obj:`float`
        Minimum ``y`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_x : :obj:`float`
        Maximum ``x`` to consider.
    max_y : :obj:`float`
        Maximum ``y`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DensityProfile"


@script_interface_register
class DipoleMoment(Observable):

    """Calculates the dipole moment for particles with given ids.

    Output format: :math:`\\left(\\sum_i q_i r^x_i, \\sum_i q_i r^y_i, \\sum_i q_i r^z_i\\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DipoleMoment"


@script_interface_register
class FluxDensityProfile(ProfileObservable):

    """Calculates the particle flux density for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider.
    min_y : :obj:`float`
        Minimum ``y`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_x : :obj:`float`
        Maximum ``x`` to consider.
    max_y : :obj:`float`
        Maximum ``y`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the x, y and z
        components of the flux density.

    """
    _so_name = "Observables::FluxDensityProfile"


@script_interface_register
class ForceDensityProfile(ProfileObservable):

    """Calculates the force density profile for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider.
    min_y : :obj:`float`
        Minimum ``y`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_x : :obj:`float`
        Maximum ``x`` to consider.
    max_y : :obj:`float`
        Maximum ``y`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the x, y and z
        components of the force.

    """
    _so_name = "Observables::ForceDensityProfile"


@script_interface_register
class LBVelocityProfile(ProfileObservable):

    """Calculates the LB fluid velocity profile.

    This observable samples the fluid in on a regular grid defined by the variables
    ``sampling_*``. Note that a small delta leads to a large number of sample
    points and carries a performance cost.

    .. WARNING::
        In case of the CPU version of the LB fluid implementation, this observable
        currently only works for a single core.

    Parameters
    ----------
    n_x_bins : :obj:`int`
        Number of bins in ``x`` direction.
    n_y_bins : :obj:`int`
        Number of bins in ``y`` direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_x : :obj:`float`
        Minimum ``x`` to consider.
    min_y : :obj:`float`
        Minimum ``y`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_x : :obj:`float`
        Maximum ``x`` to consider.
    max_y : :obj:`float`
        Maximum ``y`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.
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

    Returns
    -------
    (``n_x_bins``, ``n_y_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the x, y and z
        components of the LB velocity.

    """
    _so_name = "Observables::LBVelocityProfile"


@script_interface_register
class LBFluidStress(Observable):

    """Calculates the average stress of the LB fluid for all nodes.

    Parameters
    ----------
    None

    Returns
    -------
    (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::LBFluidStress"


@script_interface_register
class MagneticDipoleMoment(Observable):

    """Calculates the magnetic dipole moment for particles with given ids.

    Output format: :math:`\\left(\\sum_i \\mu^x_i, \\sum_i \\mu^y_i, \\sum_i \\mu^z_i\\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::MagneticDipoleMoment"


@script_interface_register
class ParticleAngularVelocities(Observable):

    """Calculates the angular velocity (omega) in the spaced-fixed frame of reference

    Output format: :math:`(\\omega^x_1,\\ \\omega^y_1,\\ \\omega^z_1),\\ (\\omega^x_2,\\ \\omega^y_2,\\ \\omega^z_2), \\dots,\\ (\\omega^x_n,\\ \\omega^y_n,\\ \\omega^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleAngularVelocities"


@script_interface_register
class ParticleBodyAngularVelocities(Observable):

    """Calculates the angular velocity (omega) in the particles'  body-fixed frame of reference.

    For each particle, the body-fixed frame of reference is obtained from the particle's
    orientation stored in the quaternions.

    Output format: :math:`(\\omega^x_1,\\ \\omega^y_1,\\ \\omega^z_1),\\ (\\omega^x_2,\\ \\omega^y_2,\\ \\omega^z_2), \\dots,\\ (\\omega^x_n,\\ \\omega^y_n,\\ \\omega^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleBodyAngularVelocities"


@script_interface_register
class ParticleBodyVelocities(Observable):

    """Calculates the particle velocity in the particles'  body-fixed frame of reference.

    For each particle, the body-fixed frame of reference is obtained from the particle's
    orientation stored in the quaternions.

    Output format: :math:`(v^x_1,\\ v^y_1,\\ v^z_1),\\ (v^x_2,\\ v^y_2,\\ v^z_2),\\ \\dots,\\ (v^x_n,\\ v^y_n,\\ v^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleBodyVelocities"


@script_interface_register
class ParticleForces(Observable):

    """Calculates the particle forces for particles with given ids.

    Output format: :math:`(f^x_1,\\ f^y_1,\\ f^z_1),\\ (f^x_2,\\ f^y_2,\\ f^z_2),\\ \\dots,\\ (f^x_n,\\ f^y_n,\\ f^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleForces"


@script_interface_register
class ParticlePositions(Observable):

    """Calculates the particle positions for particles with given ids.

    Output format: :math:`(x_1,\\ y_1,\\ z_1),\\ (x_2,\\ y_2,\\ z_2),\\ \\dots,\\ (x_n,\\ y_n,\\ z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticlePositions"


@script_interface_register
class ParticleVelocities(Observable):

    """Calculates the particle velocities for particles with given ids.

    Output format: :math:`(v^x_1,\\ v^y_1,\\ v^z_1),\\ (v^x_2,\\ v^y_2,\\ v^z_2),\\ \\dots,\\ (v^x_n,\\ v^y_n,\\ v^z_n)`.

    The particles are ordered according to the list of ids passed to the observable.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleVelocities"


@script_interface_register
class ParticleDistances(Observable):

    """Calculates the distances between particles with given ids along a
    polymer chain.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N - 1,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::ParticleDistances"


@script_interface_register
class TotalForce(Observable):

    """Calculates the total force on particles with given ids.

    Note that virtual sites are not included since forces on them do not enter the equation of motion directly.

    Output format: :math:`\\left(\\sum_i f^x_i, \\sum_i f^y_i, \\sum_i f^z_i\\right)`

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::TotalForce"


@script_interface_register
class BondAngles(Observable):

    """Calculates the angles between bonds of particles with given ids along a
    polymer chain.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N - 2,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::BondAngles"


@script_interface_register
class CosPersistenceAngles(Observable):

    """Calculates the cosine of mutual bond angles for chained particles with given ids.

    The *i*-th  value of the result contains the cosine of the angle between bonds that
    are separated by *i* bonds. The values are averaged over the chain.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N - 2,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::CosPersistenceAngles"


@script_interface_register
class BondDihedrals(Observable):

    """Calculates the dihedrals between particles with given ids along a
    polymer chain.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.

    Returns
    -------
    (N - 3,) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::BondDihedrals"


@script_interface_register
class StressTensor(Observable):

    """Calculates the total stress tensor. See :ref:`stress tensor`.

    Returns
    -------
    (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::StressTensor"


@script_interface_register
class DPDStress(Observable):

    """Calculates the non-equilibrium contribution of the DPD interaction
    to the stress tensor.

    Parameters
    ----------
    None

    Returns
    -------
    (3, 3) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::DPDStress"


@script_interface_register
class CylindricalDensityProfile(ProfileObservable):

    """Calculates the particle density in cylindrical coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``) :obj:`ndarray` of :obj:`float`

    """
    _so_name = "Observables::CylindricalDensityProfile"


@script_interface_register
class CylindricalFluxDensityProfile(ProfileObservable):

    """Calculates the particle flux density in cylindrical coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the radial distance,
        azimuth and axial coordinate of the particle flux density field.

    """
    _so_name = "Observables::CylindricalFluxDensityProfile"


@script_interface_register
class CylindricalLBFluxDensityProfileAtParticlePositions(ProfileObservable):

    """Calculates the LB fluid flux density at the particle positions in
    cylindrical coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the radial distance,
        azimuth and axial coordinate of the LB flux density field.

    """
    _so_name = "Observables::CylindricalLBFluxDensityProfileAtParticlePositions"


@script_interface_register
class CylindricalLBVelocityProfileAtParticlePositions(ProfileObservable):

    """Calculates the LB fluid velocity at the particle positions in
    cylindrical coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the radial distance,
        azimuth and axial coordinate of the LB velocity field.

    """
    _so_name = "Observables::CylindricalLBVelocityProfileAtParticlePositions"


@script_interface_register
class CylindricalVelocityProfile(ProfileObservable):

    """Calculates the particle velocity profile in cylindrical coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
        The ids of (existing) particles to take into account.
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the radial distance,
        azimuth and axial coordinate of the particle velocity field.

    """
    _so_name = "Observables::CylindricalVelocityProfile"


@script_interface_register
class CylindricalLBVelocityProfile(ProfileObservable):

    """Calculates the LB fluid velocity profile in cylindrical coordinates.

    This observable samples the fluid in on a regular grid defined by variable
    ``sampling_density``. Note that a small delta leads to a large number of
    sample points and carries a performance cost.

    Parameters
    ----------
    center : (3,) array_like of :obj:`float`
        Position of the center of the cylindrical coordinate system for the histogram.
    axis : (3,) array_like of :obj:`float`
        Orientation vector of the ``z``-axis of the cylindrical coordinate system for the histogram.
    n_r_bins : :obj:`int`
        Number of bins in radial direction.
    n_phi_bins : :obj:`int`
        Number of bins for the azimuthal direction.
    n_z_bins : :obj:`int`
        Number of bins in ``z`` direction.
    min_r : :obj:`float`
        Minimum ``r`` to consider.
    min_phi : :obj:`float`
        Minimum ``phi`` to consider.
    min_z : :obj:`float`
        Minimum ``z`` to consider.
    max_r : :obj:`float`
        Maximum ``r`` to consider.
    max_phi : :obj:`float`
        Maximum ``phi`` to consider.
    max_z : :obj:`float`
        Maximum ``z`` to consider.
    sampling_density : :obj:`float`
        Samples per unit volume for the LB velocity interpolation.

    Returns
    -------
    (``n_r_bins``, ``n_phi_bins``, ``n_z_bins``, 3) :obj:`ndarray` of :obj:`float`
        The fourth component contains the histogram for the radial distance,
        azimuth and axial coordinate of the LB velocity field.

    """
    _so_name = "Observables::CylindricalLBVelocityProfile"
