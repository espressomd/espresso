from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register



@script_interface_register
class Observable(ScriptInterfaceHelper):
    _so_name = "Observables::Observable"
    _so_bind_methods = ("calculate", "n_values")
    _so_creation_policy = "LOCAL"


@script_interface_register
class ComForce(Observable):
    """Calculates the total force on particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ComForce"


@script_interface_register
class ComPosition(Observable):
    """Calculates the center of mass for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ComPosition"


@script_interface_register
class ComVelocity(Observable):
    """Calculates the center of mass velocity for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ComVelocity"


@script_interface_register
class Current(Observable):
    """Calculates the electric current for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::Current"


@script_interface_register
class DensityProfile(Observable):
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

    """
    _so_name = "Observables::DensityProfile"


@script_interface_register
class DipoleMoment(Observable):
    """Calculates the dipole moment for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::DipoleMoment"


@script_interface_register
class FluxDensityProfile(Observable):
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

    """
    _so_name = "Observables::FluxDensityProfile"


@script_interface_register
class ForceDensityProfile(Observable):
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

    """
    _so_name = "Observables::ForceDensityProfile"


@script_interface_register
class LBVelocityProfile(Observable):
    """Calculates the LB fluid velocity profile.

    This observable samples the fluid in on a regular grid defined by the variables
    ``sampling*``. Note that a small delta leads to a large number of sample
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
        Wether or not to allow bins that will not be sampled at all.

    """
    _so_name = "Observables::LBVelocityProfile"


@script_interface_register
class MagneticDipoleMoment(Observable):
    """Calculates the magnetic dipole moment for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::MagneticDipoleMoment"


@script_interface_register
class ParticleAngularVelocities(Observable):
    _so_name = "Observables::ParticleAngularVelocities"
    
    """Calculates the angular velocity (omega) in the spaced-fixed frame of reference
   
    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """


@script_interface_register
class ParticleBodyAngularVelocities(Observable):
    _so_name = "Observables::ParticleBodyAngularVelocities"
    """Calculates the angular velocity (omega) in the particles'  body-fixed frame of reference.
   
   For each particle, the body-fixed frame of reference is obtained from the particle's
   orientation stored in the quaternions.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """


@script_interface_register
class ParticleBodyVelocities(Observable):
    """Calculates the particle velocity in the particles'  body-fixed frame of reference.
   
    For each particle, the body-fixed frame of reference is obtained from the particle's
    orientation stored in the quaternions.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ParticleBodyVelocities"


@script_interface_register
class ParticleForces(Observable):
    """Calculates the particle forces for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ParticleForces"


@script_interface_register
class ParticlePositions(Observable):
    """Calculates the particle positions for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ParticlePositions"


@script_interface_register
class ParticleVelocities(Observable):
    """Calculates the particle velocities for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ParticleVelocities"


@script_interface_register
class StressTensor(Observable):
    _so_name = "Observables::StressTensor"

    """Calculates the total stress tensor. See :ref:`stress tensor`)

    """



@script_interface_register
class CylindricalDensityProfile(Observable):
    """Calculates the particle density in polar coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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

    """
    _so_name = "Observables::CylindricalDensityProfile"


@script_interface_register
class CylindricalFluxDensityProfile(Observable):
    """Calculates the particle flux density in polar coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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

    """
    _so_name = "Observables::CylindricalFluxDensityProfile"


@script_interface_register
class CylindricalLBFluxDensityProfileAtParticlePositions(Observable):
    """Calculates the LB fluid flux density at the particle positions in polar coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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

    """
    _so_name = "Observables::CylindricalLBFluxDensityProfileAtParticlePositions"


@script_interface_register
class CylindricalLBVelocityProfileAtParticlePositions(Observable):
    """Calculates the LB fluid velocity at the particle positions in polar coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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

    """
    _so_name = "Observables::CylindricalLBVelocityProfileAtParticlePositions"


@script_interface_register
class CylindricalVelocityProfile(Observable):
    """Calculates the particle velocity profile in polar coordinates.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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

    """
    _so_name = "Observables::CylindricalVelocityProfile"


@script_interface_register
class CylindricalLBVelocityProfile(Observable):
    """Calculates the LB fluid velocity profile in polar coordinates.

    This observable samples the fluid in on a regular grid defined by the variables
    ``sampling*``. Note that a small delta leads to a large number of sample
    points and carries a performance cost.

    Parameters
    ----------
    center : array_like of :obj:`float`
             Position of the center of the polar coordinate system for the histogram.
    axis : :obj:`str` (``x``, ``y``, or ``z``)
           Orientation of the ``z``-axis of the polar coordinate system for the histogram.
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
        Wether or not to allow bins that will not be sampled at all.

    """
    _so_name = "Observables::CylindricalLBVelocityProfile"
