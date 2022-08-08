#
# Copyright (C) 2013-2022 The ESPResSo project
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

from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class GalileiTransform(ScriptInterfaceHelper):
    """
    Collective operations on particles.

    Methods
    -------
    kill_particle_motion()
        Stop the motion of all particles.

        Parameters
        ----------
        rotation : :obj:`bool`, optional
            Whether or not to stop the angular momentum too. Defaults to false.

    kill_particle_forces()
        Set the forces on the particles to zero.

        Parameters
        ----------
        torque : :obj:`bool`, optional
            Whether or not to set the torques to zero too. Defaults to false.

    system_CMS()
        Calculate the center of mass of the system. Assumes equal unit mass
        if the ``MASS`` feature is not compiled in.

        Returns
        -------
        (3,) array_like of :obj:`float`
            The center of mass position.

    system_CMS_velocity()
        Calculate the center of mass velocity of the system. Assumes equal unit
        mass if the ``MASS`` feature is not compiled in.

        Returns
        -------
        (3,) array_like of :obj:`float`
            The of the center of mass velocity vector.

    galilei_transform()
        Remove the center of mass velocity of the system. Assumes equal unit
        mass if the ``MASS`` feature is not compiled in. This is often used
        when switching from Langevin Dynamics to lattice-Boltzmann. This is
        due to the random nature of LD that yield a non-zero net system
        momentum at any given time.

    """
    _so_name = "Galilei::Galilei"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "kill_particle_motion",
        "kill_particle_forces",
        "system_CMS",
        "system_CMS_velocity",
        "galilei_transform",
    )
