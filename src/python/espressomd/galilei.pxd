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
cdef extern from "galilei.hpp":

    void local_kill_particle_motion(int rotation)
    void local_kill_particle_forces(int torque)
    void local_system_CMS(double * sdata)
    void local_system_CMS_velocity(double * sdata)
    void local_galilei_transform(double * sdata)

    ctypedef struct galilei_struct "GalileiStruct":
        double cms[3]
        double cms_vel[3]


    extern galilei_struct gal

cdef extern from "communication.hpp":

    void mpi_kill_particle_motion(int rotation)
    void mpi_kill_particle_forces(int torque)
    void mpi_system_CMS()
    void mpi_system_CMS_velocity()
    void mpi_galilei_transform() 
