from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class AutoUpdateObservables(ScriptInterfaceHelper):
    _so_name = "Observables::AutoUpdateObservables"
    _so_creation_policy = "LOCAL"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Observable):
                observable = args[0]
            else:
                raise TypeError(
                    "Either a Observable object or key-value pairs for the parameters of a Observable object need to be passed.")
        else:
            observable = Observable(**kwargs)
        self.call_method("add", object=observable)
        return observable

    def remove(self, observable):
        self.call_method("remove", object=observable)


@script_interface_register
class Observable(ScriptInterfaceHelper):
    _so_name = "Observables::Observable"
    _so_bind_methods = ("value", "calculate", "update", "auto_write_to")
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


@script_interface_register
class ParticleBodyAngularVelocities(Observable):
    _so_name = "Observables::ParticleBodyAngularVelocities"


@script_interface_register
class ParticleBodyVelocities(Observable):
    _so_name = "Observables::ParticleBodyVelocities"


@script_interface_register
class ParticleCurrent(Observable):
    """Calculates the particle current for particles with given ids.

    Parameters
    ----------
    ids : array_like of :obj:`int`
          The ids of (existing) particles to take into account.

    """
    _so_name = "Observables::ParticleCurrent"


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


@script_interface_register
class StressTensorAcf(Observable):
    _so_name = "Observables::StressTensorAcf"


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
    """Calculates the LB particle velocity profile in polar coordinates.

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
