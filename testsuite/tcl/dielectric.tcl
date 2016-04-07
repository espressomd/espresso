# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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

require_feature "ELECTROSTATICS"
require_feature "CONSTRAINTS"
require_feature "EXTERNAL_FORCES"

puts "---------------------------------------------------------------"
puts "- Testcase dielectric_test.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"


setmd box_l 20 20 20

set dist 0.5 

proc setup_wall { } {
  global n_induced_charges
  global dist
  puts "Creating Wall 1"
  dielectric wall normal 1. 0 0 dist 0.5 res 1.2
  constraint wall normal 1. 0 0 dist -0.0 type 0
  puts "Creating Wall 2"
  dielectric wall normal -1. 0 0 dist -19 res 1.2
  constraint wall normal -1. 0 0 dist -19.5 type 0
}

proc setup_sphere { } { 
  global dist
  puts "Creating sphere"
  set r 3
  dielectric sphere center 14 10 10 radius $r res 0.5 eps 0.2 
  constraint sphere center 14 10 10 radius [ expr $r-$dist ] direction 1 type 0
}

proc setup_pore { } {
  set dist 0.5
  puts "Creating pore"
  set r1 5. 
  set r2 8.
  set rs 2.5
  set l 4.
  dielectric pore  center 10 10 10 axis 1 0 0 radii $r1 $r2 length $l res 1.9 smoothing_radius $rs eps 0.2
  constraint pore center 10 10 10 axis 1 0 0 radii [ expr $r1 + $dist ] [ expr $r2 + $dist ] length [ expr $l - $dist ] type 0 smoothing_radius [ expr $rs - $dist ]
}

proc setup_slitpore { } {
  ## This tests the general case, where we prescribe the shape of the slitpore we want to get, assuming a repulsive Lennard-Jones
  ## potential with sigma=wall_sigma. All widths arguments given to constraints are made larger accordingly and smoothing radii are adapted.
  global dist
  set pore_mouth 10
  set d $dist
  set channel_width 4.
  set pore_width 3.
  set pore_length 4. 
  set upper_smoothing_radius 2.
  set lower_smoothing_radius 2. 
  set wall_sigma 1.
  
  dielectric slitpore pore_mouth [ expr $pore_mouth - $d ] \
                      channel_width [ expr $channel_width+2*$d] \
                      pore_width [ expr $pore_width + 2*$d] \
                      pore_length [ expr $pore_length + 0*$d ] \
                      upper_smoothing_radius [ expr $upper_smoothing_radius -$d ] \
                      lower_smoothing_radius [ expr $lower_smoothing_radius + $d ] \
                      res 0.5 eps 100
  
  constraint slitpore pore_mouth [ expr $pore_mouth - $wall_sigma ] \
                      channel_width [ expr $channel_width+2*$wall_sigma] \
                      pore_width [ expr $pore_width + 2*$wall_sigma] \
                      pore_length [ expr $pore_length + 0*$wall_sigma ] \
                      upper_smoothing_radius [ expr $upper_smoothing_radius -$wall_sigma ] \
                      lower_smoothing_radius [ expr $lower_smoothing_radius + $wall_sigma ] \
                      type 0
}


#dielectric pore center 10 10 10 axis 1 0 0 radii $r1 $r2 length $l res 1.9 smoothing_radius $rs eps 0.2



proc setup_cylinder { } {
  puts "Creating cylinder"
  set r 2
  set l 3
  dielectric_cylinder center 10 10 10   axis 1 0 0 radius $r length $l res .25 covers 1
  constraint cylinder center 10 10 10   axis 1 0 0 radius [ expr $r - $dist ] length [ expr $l - $dist ] type 0 direction 1
  puts "Done. Now we have $n_induced_charges induced charges"
}

proc check_distances { } {
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  set dist 0.5
  puts "Checking all distances and normal vectors"
  for { set i 0 } { $i < $n_induced_charges } { incr i } {
    set pos [ part $i print pos ]
    set distance [ constraint mindist_position [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
    set dv [ constraint mindist_position_vec [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
    set n  [ lindex $icc_normals $i ] 
    set diff [ expr [ lindex $dv 0 ] /$dist - [ lindex $n 0 ]] 
    if { $diff > 1e-6 } {
      puts "Everything is wrong at id $i!"
      puts "pos is $pos"
      puts "dist is $distance"
      puts "dist_vec is $dv"
      puts "normal is $n"
      return 1
    }
  }
  puts "OK"
}

proc reset { } {
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  for { set i 0 } { $i < $n_induced_charges } { incr i } {
    part $i delete 
  }
  part 0 pos 1 2 3
  part 0 delete
  set n_induced_charges 0 
  unset icc_areas 
  unset icc_normals
  unset icc_epsilons
  unset icc_sigmas
  for { set i  [ constraint n_constraints ] } { $i > 0 } { set i [ expr $i -1 ] } {
    constraint delete  [ expr $i-1 ]
  }
}

##### We can enable this to make visual checks of everything!

proc write_distfile { fname } {
  set distfile [ open $fname "w" ]
  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    set pos [ part $i print pos ]
    set dv [ constraint mindist_position_vec [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
    puts $distfile "$pos $dv" 
  }
  close $distfile
}

#
setup_wall 
puts "Done. Now we have $n_induced_charges induced charges"
puts "Espresso believes we have [ setmd n_part ] particles"
check_distances
reset
#
setup_sphere
puts "Done. Now we have $n_induced_charges induced charges"
puts "Espresso believes we have [ setmd n_part ] particles"
check_distances
reset
#

setup_slitpore
puts "Done. Now we have $n_induced_charges induced charges"
puts "Espresso believes we have [ setmd n_part ] particles"
check_distances
reset


setup_pore
puts "Done. Now we have $n_induced_charges induced charges"
puts "Espresso believes we have [ setmd n_part ] particles"
check_distances
reset
#prepare_vmd_connection "test" 3000
#after 200000

