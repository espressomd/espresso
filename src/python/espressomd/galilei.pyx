#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

from __future__ import print_function, absolute_import
from . cimport galilei

cdef class GalileiTransform:

    def kill_particle_motion(self, rotation = 0):
        """ Stop motion of the particles """
        mpi_kill_particle_motion(rotation)

    def kill_particle_forces(self, torque = 0):
        """ Set the forces on the particles to zero """
        mpi_kill_particle_forces(torque)

    def system_CMS(self):
        """ Calculate the CMS of the system"""
        mpi_system_CMS()
        return gal.cms

    def system_CMS_velocity(self):
        """ Calculate the CMS velocity of the system"""
        mpi_system_CMS_velocity()
        return gal.cms_vel

    def galilei_transform(self):
        """ Remove the CMS velocity of the system """
        mpi_galilei_transform()
