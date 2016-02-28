# Copyright (C) 2011,2012,2013,2016 The ESPResSo project
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

### Measuring the force on a single sphere immersed in a fluid with
### fixed velocity boundary conditions created by two 
### walls at finite distance.
### The force is compared to th analytical result F=6 pi eta r v
### i.e. the stokes force on the particles.


# We create a box of size box_width x box_width x box_length and 
# place an object in the center. We measure the drag force 
# in z direction. We create walls in the xz and yz plane at the box 
# boundaries, where the velocity is fixed to $v. 
#
set box_width 64
set box_length 64


# The temperature is zero.
thermostat lb 0.
# Espresso complains if we don't give a skin
setmd skin 0.5

# The boundary speed 
set v 0.01

# The fluid parameters
set kinematic_visc 1.
set density 1.
set agrid 1.
set tau 0.1
# not necessary, but LB complains without
set friction 10. 

# Setting box length
setmd box_l [ expr $box_width+2*$agrid ]  [ expr $box_width+2*$agrid ] $box_length 
# MD timestep and LB timestep are identical 
setmd time_step $tau
# Invoke LB fluid
lbfluid gpu visc $kinematic_visc dens $density  agrid $agrid tau $tau friction $friction

# Four walls make an infinite square channel along z direction
lbboundary wall normal -1. 0. 0. dist [ expr -(1+$box_width) ] velocity 0.00 0 $v 
lbboundary wall normal  1. 0. 0. dist 1. velocity 0. 0. $v
lbboundary wall normal  0 -1. 0. dist [ expr -(1+$box_width) ] velocity 0.00 0 $v 
lbboundary wall normal  0  1. 0. dist 1. velocity 0. 0. $v


set radius 5.5
lbboundary sphere center [ expr 0.5*($box_width+2*$agrid) ] [ expr 0.5*($box_width+2*$agrid) ] [ expr 0.5*$box_length ]\
     radius $radius direction +1

lbfluid print vtk boundary "boundary.vtk"
for { set i 0 } { $i < 30 } { incr i } {
  # integrate a number of steps
  integrate 400
  # write flow field to disk
  lbfluid print vtk velocity [ format "fluid%04d.vtk" $i ]
  # print out force on the sphere
  puts [ lindex [  lbboundary force  4 ]  2 ]
}

set stokes_force [ expr 6*3.1415*$kinematic_visc*$radius*$v ]
puts "Stokes' Law says: f=$stokes_force"
