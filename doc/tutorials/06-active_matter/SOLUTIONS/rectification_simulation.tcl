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
#                     Active Matter: Rectification Tutorial                    #
#                                                                              #
################################################################################

source "tutorial_tests.tcl"

require_feature "ENGINE"
require_feature "CONSTRAINTS"
require_feature "LENNARD_JONES"
require_feature "ROTATION"
require_feature "MASS"

# Quaternion procedure

proc a2quat {phi theta} {
  
  set q1w   [expr cos($theta/2.0)]
  set q1x   0
  set q1y   [expr sin($theta/2.0)]
  set q1z   0

  set q2w   [expr cos($phi/2.0)]
  set q2x   0
  set q2y   0
  set q2z   [expr sin($phi/2.0)]

  set q3w   [expr ($q1w*$q2w-$q1x*$q2x-$q1y*$q2y-$q1z*$q2z)]
  set q3x   [expr ($q1w*$q2x+$q1x*$q2w-$q1y*$q2z+$q1z*$q2y)]
  set q3y   [expr ($q1w*$q2y+$q1x*$q2z+$q1y*$q2w-$q1z*$q2x)]
  set q3z   [expr ($q1w*$q2z-$q1x*$q2y+$q1y*$q2x+$q1z*$q2w)]
  return    [list $q3w $q3x $q3y $q3z]
}

################################################################################

# Read in the active velocity from the command prompt

if { [llength $argv] != 1 } {
  puts "Usage: Espresso $argv0 <vel> (0 <= vel < 10.0)"
  exit 1
}

set vel [lindex $argv 0]

################################################################################

# Setup constants

set PI [expr acos(-1.0)]
set tcl_precision 8

set outdir "RESULTS_RECTIFICATION"
file mkdir $outdir
after 250

# Setup the box (we pad the diameter to ensure that the LB boundaries
# and therefore the constraints, are away from the edge of the box)

set length 100
set diameter 20
setmd box_l $length [expr $diameter + 4] [expr $diameter + 4]

# Setup the MD parameters (Langevin -> no hydrodynamics)

set prod_steps 500
set prod_length 500
set dt 0.01
setmd periodic 1 1 1
setmd time_step $dt
setmd skin 0.1
setmd min_global_cut 0.5
thermostat langevin 1.0 1.0

################################################################################
# 
# Here we use exactly the same parameters for the geometry of the constraints
# that was used for the LB boundaries. This can be done, since the distance
# function used for the constraints is the same as the one used for the 
# LB boundaries.
#
################################################################################

# Setup cylinder

constraint cylinder \
  center [expr $length/2.0] \
         [expr ($diameter+4)/2.0] \
         [expr ($diameter+4)/2.0] \
  axis 1 0 0 \
  radius [expr $diameter/2.0] \
  length [expr $length] \
  direction -1 \
  type 1

# Setup walls

constraint wall dist 2 normal 1 0 0 type 2
constraint wall dist [expr -($length - 2)] normal -1 0 0 type 3

# Setup cone

set irad 4.0
set angle [expr $PI/4.0]
set orad [expr ($diameter - $irad)/sin($angle)]
set shift [expr 0.25*$orad*cos($angle)]

constraint hollow_cone \
  center [expr $length/2.0 - $shift] \
         [expr ($diameter+4)/2.0] \
         [expr ($diameter+4)/2.0] \
  orientation 1 0 0 \
  outer_radius $orad \
  inner_radius $irad \
  width 2.0 \
  opening_angle $angle \
  direction 1 \
  type 4

################################################################################
#
# We set up a WCA (almost-hard) interaction between the particles and the 
# the confining geometry. We do not have particle-particle interactions, which
# are not necessary to observe rectification. 
#
################################################################################

set sig 0.5
set cut [expr 1.12246*$sig]
set eps 1.0
set shift [expr 0.25]

inter 0 1 lennard-jones $eps $sig $cut $shift 0
inter 0 2 lennard-jones $eps $sig $cut $shift 0
inter 0 3 lennard-jones $eps $sig $cut $shift 0
inter 0 4 lennard-jones $eps $sig $cut $shift 0

################################################################################
#
# Setup the particles. We put them all in two points one in each chamber
# and give them random directions. This speeds up the equilibration, since
# putting them all in a single chamber, would make it take a long time to
# observe the effect of rectification.
#
################################################################################

set cntr 0
set npart 500
while {$cntr < $npart} {

  if {$cntr % 2 == 0} {
    set x [expr 0.25*$length]
  } else {
    set x [expr 0.75*$length]
  }

  set y [expr ($diameter+4)/2.0]
  set z [expr ($diameter+4)/2.0]

  set theta [expr 2*[t_random]*$PI]
  set phi   [expr 2*[t_random]*$PI]
  set quats [a2quat $theta $phi]

  set idx [setmd n_part]
  part $idx pos $x $y $z type 0 swimming v_swim $vel \
  quat [lindex $quats 0] [lindex $quats 1] [lindex $quats 2] [lindex $quats 3]

  incr cntr
}

################################################################################

# Equilibrate

integrate [expr 25*$prod_length]

# Output the CMS coordinates

set outfile [open "$outdir/CMS\_$vel.dat" "w"]
puts $outfile "####################################################"
puts $outfile "#        time    CMS x coord    average CMS        #"
puts $outfile "####################################################"
flush $outfile

# Production run

set dev_sum 0.0
set dev_av 0.0
for { set i 0 } { $i < $prod_steps} { incr i } {

  # We output the coordinate of the center of mass in 
  # the direction of the long axis, here we consider 
  # the deviation from the center

  set dev [expr [lindex [system_CMS] 0] - 0.5*$length];

  if { $i > 0 } {
    set dev_sum [expr $dev_sum + $dev]
    set dev_av [expr $dev_sum/$i]
  }

  puts $outfile "[setmd time] $dev $dev_av"
  flush $outfile

  integrate $prod_length
}

close $outfile

# Output the final configuration

writevtk "$outdir/points\_$vel.vtk" 0

################################################################################

exit 0

################################################################################
