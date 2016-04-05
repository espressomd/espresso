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
if {1} {
set box_l_x 20
set box_l_y 4
set box_l_z 30


setmd box_l $box_l_x $box_l_y $box_l_z
set vmd_output "yes"

set pi 3.14144
set res 1.0

set pore_mouth 20.
set channel_width 5.
set box_l_x 30.
set pore_length 10
set pore_width 0.8
set upper_smoothing_radius 1.000000
set lower_smoothing_radius 1.0



set d 1.
#inter 0 0 lennard-jones 1. 1.0 1.12245 0.5 0.0
inter 1 2 lennard-jones 1. 1.0 1.12245 0.5 0.0
#inter 1 1 lennard-jones 1. 1.0 1.12245 0.5 0.0


constraint slitpore pore_mouth [ expr $pore_mouth - $d ] channel_width [ expr $channel_width+2*$d] pore_width [ expr $pore_width + 2*$d] pore_length [ expr $pore_length + 0*$d ] upper_smoothing_radius [ expr $upper_smoothing_radius -$d ] lower_smoothing_radius [ expr $lower_smoothing_radius + $d ] type 1
#constraint wall normal 1 0 0 dist 5 type 0
#constraint wall normal -1 0 0 dist -25 type 0
#constraint slitpore pore_mouth $pore_mouth channel_width $channel_width pore_length $pore_length pore_width $pore_width upper_smoothing_radius $upper_smoothing_radius lower_smoothing_radius $lower_smoothing_radius
set d 0.5
set n_induced_charges 0
dielectric slitpore pore_mouth [ expr $pore_mouth - $d ] channel_width [ expr $channel_width+2*$d] pore_width [ expr $pore_width + 2*$d] pore_length [ expr $pore_length + 0*$d ] upper_smoothing_radius [ expr $upper_smoothing_radius -$d ] lower_smoothing_radius [ expr $lower_smoothing_radius + $d ] res 0.5

global n_induced_charges
set partcounter $n_induced_charges

for { set  i 0 } { $i < 200 } { incr i } {
  set ok 0 
  while { !$ok } {
    set x [ expr [ t_random ] * $box_l_x ]
    set y [ expr [ t_random ] * $box_l_y ]
    set z [ expr [ t_random ] * $box_l_z ]
  
    if { [ constraint mindist_position $x $y $z ] > 2.3  } { 
      set ok 1
    }
  }
  part $partcounter pos $x $y $z type  2 
  incr partcounter
}

setmd time_step 0.01
thermostat langevin 1. 1. 
setmd skin 0.5
integrate 0


if { $vmd_output == "yes" } {
  prepare_vmd_connection "test" 10000 
  after 1000
  imd positions
}

for { set i 0 } { $i < 1000 } { incr i } {
  integrate 100
  imd positions
}



#set distfile [ open "dist.dat" "w" ]
#for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
#  set pos [ part $i print pos ]
#  set dv [ constraint mindist_position [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
#  puts $distfile "$pos $dv" 
#}
#close $distfile



set ofile [ open "test.vtf" "w" ]
writevsf $ofile
writevcf $ofile
} 
exec killall vmd

