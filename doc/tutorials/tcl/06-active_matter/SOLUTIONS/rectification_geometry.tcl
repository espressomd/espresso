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
#                   Active Matter: Rectification System Setup                  #
#                                                                              #
################################################################################

source "tutorial_tests.tcl"

require_feature "LB_GPU"
require_feature "LB_BOUNDARIES_GPU"

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

# Setup the MD parameters 

set dt 0.01
setmd periodic 1 1 1
setmd time_step $dt
setmd skin 0.1
setmd min_global_cut 0.5

# Setup LB parameters (these are irrelevant here) and fluid

set agrid 1
set vskin 0.1
set frict 20.0
set visco 1.0
set densi 1.0

lbfluid gpu agrid $agrid dens $densi visc $visco tau $dt friction $frict

################################################################################
# 
# Now we set up the three LB boundaries that form the rectifying geometry.
# The cylinder boundary/constraint is actually already capped, but we put
# in two planes for safety's sake. If you want to create an cylinder of
# 'infinite length' using the periodic boundaries, then the cylinder must
# extend over the boundary.
#
################################################################################

# Setup cylinder

lbboundary cylinder \
  center [expr $length/2.0] \
         [expr ($diameter+4)/2.0] \
         [expr ($diameter+4)/2.0] \
  axis 1 0 0 \
  radius [expr $diameter/2.0] \
  length [expr $length] \
  direction -1

# Setup walls

lbboundary wall dist 2 normal 1 0 0
lbboundary wall dist [expr -($length - 2)] normal -1 0 0 

# Setup cone

set irad 4.0
set angle [expr $PI/4.0]
set orad [expr ($diameter - $irad)/sin($angle)]
set shift [expr 0.25*$orad*cos($angle)]

lbboundary hollow_cone \
  center [expr $length/2.0 - $shift] \
         [expr ($diameter+4)/2.0] \
         [expr ($diameter+4)/2.0] \
  orientation 1 0 0 \
  outer_radius $orad \
  inner_radius $irad \
  width 2.0 \
  opening_angle $angle \
  direction 1

################################################################################

# Output the geometry

lbfluid print vtk boundary "$outdir/boundary.vtk"

################################################################################

exit 0

################################################################################
