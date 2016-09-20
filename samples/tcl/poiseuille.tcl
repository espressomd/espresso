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

setmd box_l 32 32 32
set dt 0.01
setmd time_step $dt

#note that the skin below does nothing as there are no particles and hence no Verlet lists which require a skin
#the variable must be set however as otherwise Espresso throws an error
setmd skin 0.3

#The command below create two plane walls effectively at y=1 and y=31.
#This will result in boundary nodes at y=0.5 and y=31.5.
lbboundary wall dist 1.0 normal 0.0 1.0 0.0
lbboundary wall dist -31.0 normal 0.0 -1.0 0.0


# set the fluid  
#note that the time step only needs to be larger than the md time step 
lbfluid agrid 1.0 dens 1.0 visc 3.0 tau $dt ext_force 0.001 0.0 0.0 friction 20.0
#note that adding the thermostat slows the simulation significantly and thermostat lb 0.0 actually does nothing
#more appropriate would be thermostat off but thermostat lb 0.0 is left in as an example of how to set the thermostat
thermostat lb 0.0

#here we print out the boundaries once at the start and the fluid over time in vtk format, ideal for Paraview
lbfluid print vtk boundary lb_fluid_boundary.vtk
for {set i 0} {$i<100} {incr i} {
 lbfluid print vtk velocity lb_fluid_vel_frame_$i.vtk
 #note that the number passed to integrate is the number of MD time steps to integrate, though in this case it is the same as the number of LB time steps
 integrate 100
}
