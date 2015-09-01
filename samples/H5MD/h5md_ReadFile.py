#
# Copyright (C) 2013,2014 The ESPResSo project
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

from __future__ import print_function
import espressomd._system as es
import espressomd
import h5md
import numpy

# H5MD: Create/Read File dataset
#############################################################
system = espressomd.System()
h5_filename="h5md_WriteReadFile.h5"
h5=h5md.h5md("h5md_WriteReadFile.h5",system)

n_part=10
n_time=10
int_steps=10

#h5.h5_dataset_size(h5.h5_file['particles/position/step'])
h5.h5_read_particles.pos(n_time-1)
#h5.h5_read_particles.image(n_time-1)
h5.h5_read_particles.v(n_time-1)
h5.h5_read_particles.f(n_time-1)
h5.h5_read_particles.type()
h5.h5_read_particles.mass()
h5.h5_read_particles.id()
#h5.h5_read_particles.box(n_time-1)
#h5.h5_read_observable(n_time-1,[2,3,4],"energies","Obs1")

for i in range(n_part):
	print(system.part[i].pos)
















# #H5MD: Write VMD parameters
# h5.h5_write_vmd_parameters()
# h5.h5_write_vmd_parameters_extra.chain(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.name(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.resid([1,2,3,4,5])
# h5.h5_write_vmd_parameters_extra.resname(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.segid(["A", "B", "C"])
# h5.h5_write_vmd_parameters_extra.type(["A", "B", "C"])


# h5.h5_read_vmd_parameters("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.chain("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.name("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.resid("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.resname("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.segid("A")#xxxxxxxxxxxxxxxx
# h5.h5_read_vmd_parameters_extra.type("A")#xxxxxxxxxxxxxxxx


#TODO def MASS  ,  remove es aus constructor  ,  remove reload file 
#TODO enegrie string to values
#TODO image
#attributes
#box_l nur [10,0,0]
