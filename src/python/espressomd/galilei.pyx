#
# Copyright (C) 2013-2019 The ESPResSo project
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

from . cimport galilei
from .utils cimport make_array_locked

cdef class GalileiTransform:

    def kill_particle_motion(self, rotation=False):
        """
        Stop the motion of the particles.

        Parameters
        ----------
        rotation : :obj:`bool`, optional
                   Whether or not to kill the rotations too.

        """
        mpi_kill_particle_motion(int(rotation))

    def kill_particle_forces(self, torque=False):
        """
        Set the forces on the particles to zero.

        Parameters
        ----------
        torque : :obj:`bool`, optional
                 Whether or not to kill the torques on all particles too.

        """
        mpi_kill_particle_forces(int(torque))

    def system_CMS(self):
        """
        Calculate the center of mass of the system. Assumes equal unit mass if the mass feature is not used.

        Returns
        -------
        cms : (3,) array_like of :obj:`float`
              The of the center of mass position vector.

        """
        return make_array_locked(mpi_system_CMS())

    def system_CMS_velocity(self):
        """
        Calculate the center of mass velocity of the system. Assumes equal unit
        mass if the mass feature is not used.

        Returns
        -------
        cms_vel : (3,) array_like of :obj:`float`
                  The of the center of mass velocity vector.

        """

        return make_array_locked(mpi_system_CMS_velocity())

    def galilei_transform(self):
        """
        Remove the center of mass velocity of the system. Assumes equal unit
        mass if the mass feature is not used. This is often used when switching
        from Langevin Dynamics to lattice-Boltzmann. This is due to the random
        nature of LD that yield a non-zero net system momentum at any given
        time.

        """
        mpi_galilei_transform()
