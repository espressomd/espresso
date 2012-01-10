#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
exec $ESPRESSO_SOURCE/Espresso $0 $*
#
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#############################################################
#                                                           #
#  Functions (TCL procedures used in lj_tutorial )          #
#                                                           # 
#############################################################
# Velocity  Rescaling
proc rescale_velocities { target_temperature particle_number } {
        set energies [analyze energy]
        set kinetic [lindex $energies 1 1]
        set factor [expr sqrt(0.5*$target_temperature*(3.0*$particle_number-3.0)/$kinetic)]
        for {set i 0} {$i<$particle_number} {incr i} {
                set vel [part $i print v]
                part $i v [expr [lindex $vel 0]*$factor] [expr [lindex $vel 1]*$factor] [expr [lindex $vel 2]*$factor]
        }
}
proc save_sim {cfile parinfo range } {
# write all available sim information to channel cfile
# in block format
 blockfile $cfile write variable all
 blockfile $cfile write tclvariable all
 blockfile $cfile write particles $parinfo $range
 blockfile $cfile write interactions
 blockfile $cfile write bonds
 blockfile $cfile write random
 blockfile $cfile write seed
 blockfile $cfile write bitrandom
 blockfile $cfile write bitseed
}
