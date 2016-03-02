# Copyright (C) 2014,2015,2016 The ESPResSo project
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
file delete -force "output"
file mkdir output

# what input files to read
set fileNodes  "input/cell_nodes.dat"
set fileTriangles "input/cell_triangles.dat"

# integrator settings for the simulation 
setmd time_step 0.1    
setmd skin 0.4
thermostat off

# rectangular channel box geometry
set boxX 50
set boxY 22
set boxZ 20
setmd box_l $boxX $boxY $boxZ

# initialization of the object-in-fluid mechanisms
oif_init

# creating templates
oif_create_template template-id 0 nodes-file $fileNodes triangles-file $fileTriangles stretch 3.0 3.0 3.0 ks 0.07 kb 0.01 kal 0.01 kag 0.01 kv 10.0

set pi 3.14159265359

# adding cells 
oif_add_object object-id 0 template-id 0 origin 5 15 5 rotate 0 0 [expr $pi/2] part-type 0 mass 1
oif_add_object object-id 1 template-id 0 origin 5 5 15 rotate 0 0 0 part-type 1 mass 1

# cell-cell interactions 
inter 0 1 soft-sphere 0.005 2.0 0.3 0.0

# cell-wall interactions 
inter 0 10 soft-sphere 0.0001 1.2 0.1 0.0

# set up fluid
lbfluid grid 1 dens 1.0 visc 1.5 tau 0.1 friction 0.5

source boundaries.tcl

# setting the constant velocity
# of the fluid on the left side of the md_box
lbboundary rhomboid velocity 0.005 0 0 corner 0 1 1 a 1 1 1 b 0 [expr $boxY-1] 1 c 0 1 [expr $boxZ-1] direction 1 

# main iteration loop
set steps 200
set counter 0
while { $counter<200} {

    set cycle [expr $counter*$steps]
    puts "cycle $cycle" 
    lbfluid print vtk velocity "output/fluid$cycle.vtk" 
    oif_object_output object-id 0 vtk-pos "output/cell0_$cycle.vtk"
    oif_object_output object-id 1 vtk-pos "output/cell1_$cycle.vtk"
    integrate $steps
    incr counter
}
