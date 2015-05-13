# Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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
source "tests_common.tcl"



proc verify { i v o } {
  # Verify that the ith component of o is zero if v is zero or non-zero else.
  
  set x [lindex $o $i]
  if {$v ==1 } {
    if { $x == 0 } {
      error "The $i'th component of the angular velocity was zero even though it should rotate"
    }
  }
  if {$v ==0 } {
    if { $x != 0 } {
      error "The $i'th component of the angular velocity was non-zero even though it should be blocked"
    }
  }
}

require_feature "ROTATION"
require_feature "ROTATION_PER_PARTICLE"

puts "----------------------------------------------"
puts "- Testcase rotation_per_particle.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

set epsilon 5e-4
thermostat langevin 1 1

setmd time_step 0.01
setmd skin 0


part 0 pos 0 0 0

# Check that rotation is on by default and that the rotation al degrees are thermalized
set rot [part 0 print rotation]
if {$rot != "14"} {
 error "Rotation not on by default!"
}

integrate 100
set quat [part 0 print quat]
set zero "1. 0. 0. 0."
if {[veclen [vecsub $quat $zero]] < $epsilon} {
 error "Rotational degrees not thermalized when rotation is on"
}

# Check that the rotational properties don't change any more, when rotation is off
set quat [part 0 print quat]
set omega [part 0 print omega_lab]
part 0 rotation 0
integrate 1000
set quatN [part 0 print quat]
set omegaN [part 0 print omega_lab]

if {[veclen [vecsub $quat $quatN]] >$epsilon} {
 error "Quaternions changed even when rotation is off"
}

if {[veclen [vecsub $omega $omegaN]] > $epsilon} {
 error "omega changed even when rotation is off"
}



# Test blocking individual axes

# Flags to activate rotation around axes

foreach x {0 1} {
  foreach y {0 1} {
    foreach z {0 1} {
      set rot 0
      if { $x==1} {
        set rot [expr $rot +2]
      }	
      if { $y==1} {
        set rot [expr $rot +4]
      }	
      if { $z==1} {
        set rot [expr $rot +8]
      }	
      part 0 rotation $rot quat 1 0 0 0 omega_body 0 0 0 
      #ext_torque -5 2 3
      if { [part 0 print rotation] != $rot } {
        error "Rotation particle property did not get set to correct value"
      }

      integrate 100
      set o [part 0 print omega_body]
      puts "$x $y $z $rot $o"
      verify 0 $x $o
      verify 1 $y $o
      verify 2 $z $o
    }
  }
}

exit 0
