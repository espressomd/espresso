#
# Copyright (C) 2023 The ESPResSo project
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

import numpy as np
from . import code_features
from .math import calc_quaternions_from_angles
if code_features.has_features("VIRTUAL_SITES_RELATIVE"):
    from .virtual_sites import VirtualSitesRelative


def add_dipole_particle(system, particle, dipole_length,
                        dipole_particle_type, mode="pusher"):
    """
    Add a virtual particle pointing opposite of the swimmer particle
    either in front of (mode="puller") or behind (mode="pusher") the swimmer.
    The virtual particle is set up to exert the swimming force of the swimmer
    particle on a fluid in the opposite direction.

    Parameters
    ----------
    system : :obj:`~espressomd.system.System`
        The system that the particle belongs to.
    particle : :obj:`~espressomd.particle_data.ParticleHandle`
        The swimmer particle that will be paired with the virtual dipole particle.
        Must have the ``swimming`` attribute already set up to get the correct
        value for ``f_swim``.
    dipole_length : :obj:`float`
        The distance between the swimmer and the virtual dipole particle.
    dipole_particle_type : :obj:`int`
        The type of the virtual dipole particle.
    mode : :obj:`str`
        Allowed values: "pusher" (default) and "puller".
        Determines whether the virtual dipole particle
        will be placed in front of or behind the swimmer.

    Returns
    -------
    dipole_particle : :obj:`~espressomd.particle_data.ParticleHandle`
        The newly created particle.

    """
    if mode == "pusher":
        dip_sign = -1
    elif mode == "puller":
        dip_sign = 1
    else:
        raise ValueError(f"'mode' must be 'pusher' or 'puller', got '{mode}'")

    if dipole_length < 0:
        raise ValueError("'dipole_length' must be >= 0.")

    code_features.assert_features(["ENGINE", "VIRTUAL_SITES_RELATIVE"])
    if not isinstance(system.virtual_sites, VirtualSitesRelative):
        raise RuntimeError(
            "system.virtual_sites must be espressomd.virtual_sites.VirtualSitesRelative.")
    if not system.virtual_sites.have_quaternion:
        raise RuntimeError(
            "system.virtual_sites must have quaternion option turned on ('have_quaternion = True').")

    p = system.part.add(
        pos=particle.pos + dip_sign * dipole_length * particle.director,
        virtual=True, type=dipole_particle_type, swimming={
            "f_swim": particle.swimming["f_swim"],
            "is_engine_force_on_fluid": True,
        })
    p.vs_auto_relate_to(particle)
    p.vs_quat = calc_quaternions_from_angles(np.pi, 0.)
    return p
