################################################################################
#                                                                              #
# Copyright (C) 2010,2011,2012,2013,2014, 2015,2016 The ESPResSo project            #
#                                                                              #
# This file is part of ESPResSo.                                               #
#                                                                              #
# ESPResSo is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# ESPResSo is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.        #
#                                                                              #
################################################################################
#                                                                              #
#                  Active Matter: Swimmer Flow Field Tutorial                  #
#                                                                              #
################################################################################

source "tutorial_tests.tcl"

require_cudadevice
require_feature "ENGINE"
require_feature "LB_GPU"
require_feature "MASS"
require_feature "ROTATION"
require_feature "ROTATIONAL_INERTIA"

# Read in the hydrodynamic type (pusher/puller) and position

if { [llength $argv] != 2 } {
  puts "Usage: Espresso $argv0 <type> <pos>"
  exit 1
}

set type   [lindex $argv 0]
set pos    [lindex $argv 1]

################################################################################

set PI [expr acos(-1.0)]
set tcl_precision 8

set outdir "./RESULTS_FLOW_FIELD/T\_$type\_P\_$pos/"
file mkdir $outdir
after 250

# System parameters

set length 25.0
setmd box_l $length $length $length
setmd periodic 1 1 1

set prod_steps 1000
set prod_length 50
set dt 0.01
setmd time_step $dt
setmd skin 0.1
setmd min_global_cut 1.0

################################################################################

# Set the position of the particle

set x0pnt [expr 0.5*$length]
set y0pnt [expr 0.5*$length]
set z0pnt [expr 0.5*$length]

set x 0.0
set y 0.0
set z $pos

set x0 [expr $x0pnt + $x]
set y0 [expr $y0pnt + $y]
set z0 [expr $z0pnt + $z]

# Sphere size, mass, and moment of inertia, dipole force

set sph_size      0.5
set sph_mass      4.8
set Ixyz          4.8
set force         0.1 

# Setup the particle particle

set cent [setmd n_part]
part $cent pos $x0 $y0 $z0 type 0 mass $sph_mass rinertia $Ixyz $Ixyz $Ixyz \
     swimming f_swim $force $type dipole_length [expr $sph_size + 0.5]

################################################################################

# Setup the fluid (quiescent)

set agrid 1
set vskin 0.1
set frict 20.0
set visco 1.0
set densi 1.0
set temp 0.0

lbfluid gpu agrid $agrid dens $densi visc $visco \
            tau $dt friction $frict couple 3pt

thermostat lb $temp

################################################################################

# Output the coordinates

set outfile [open "$outdir/trajectory.dat" "w"]
puts $outfile "####################################################"
puts $outfile "#        time        position       velocity       #"
puts $outfile "####################################################"
flush $outfile

# Production run

for { set k 0 } { $k <= $prod_steps } { incr k } {

  # Output quantities

  puts $outfile "[setmd time] [part 0 print pos v]"
  flush $outfile

  # Output 50 simulations

  if { [expr $k % ($prod_steps/50)] == 0 } {
    set num [expr $k/($prod_steps/50)]
    lbfluid print vtk velocity "$outdir/lb_velocity\_$num.vtk" 
    writevtk "$outdir/position\_$num.vtk" 0
  }

  # Integrate

  integrate $prod_length
}

close $outfile

################################################################################

exit 0

################################################################################
