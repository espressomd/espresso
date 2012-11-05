#
# Copyright (C) 2010,2012 The ESPResSo project
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
#############################################################
#                                                           #
#  Sample System : 1:1 Electrolyte in the Restricted        #
#  Primitive Model confined to a slit pore                  #
#                                                           # 
#############################################################

# set the properties of the system. The number of particles 
# must be divisible by 4 to ensure a charge neutral system

set n_part 200
set density 0.7

# set the box length in every direction. Remember we are 
# interested in a slit-like geometry, which is small in
# one direction and periodic in the other. Note we also 
# compensate for the 'size' of the particles and the hard 
# walls in setting up the slit dimensions.

set box_l [expr pow($n_part/$density,1./3.)]
set box_lxy [expr sqrt(2)*$box_l]
set box_lz [expr 0.5*$box_l]

# we pass the box dimension to ESPResSo
setmd box_l $box_lxy $box_lxy [expr $box_lz + 1.0]
setmd periodic 1 1 0

# next we set up two 'hard' walls. A wall is defined by
# its normal and its distance to the origin
constraint wall normal 0 0 1 dist 0 type 2
constraint wall normal 0 0 -1 dist [expr -$box_lz - 1.0] type 2

# we set up two charged plates with a charge density sigma
# at the location of the walls. The charge on the plates 
# negative
set sigma [expr -0.25*$n_part/($box_lxy*$box_lxy)]
constraint plate height 0 sigma $sigma
constraint plate height [expr $box_lz + 1.0] sigma $sigma

# set the particles at random positions in the box
# note that we compensate for the wall positions in
# the z direction. Note that we first set up the
# same number of positive and negative particles. 
# Next we set up all the positive charges that are 
# required to make the system charge neutral

t_random seed 12345
set q 1
set type 0
for {set i 0} { $i < $n_part/2 } {incr i} {
  set posx [expr $box_lxy*[t_random]]
  set posy [expr $box_lxy*[t_random]]
  set posz [expr ($box_lz-1.0)*[t_random] + 1.0]
  set q [expr -$q]
  set type [expr 1-$type]
  part $i pos $posx $posy $posz q $q type $type
}

set q 1
set type 0
for {set i 0} { $i < $n_part/2 } {incr i} {
  set posx [expr $box_lxy*[t_random]]
  set posy [expr $box_lxy*[t_random]]
  set posz [expr ($box_lz-1.0)*[t_random] + 1.0]
  part [expr $i + $n_part/2] pos $posx $posy $posz q $q type $type
}

# set the properties of the integrator, we are
# interested in a proper NVT ensemble, so we use
# the Langevin thermostat

setmd time_step 0.01
setmd skin 0.3
set temp 1
set gamma 1
thermostat langevin $temp $gamma

# set the 'hard' core of the particles via a 
# Lennard-Jones potential, in this case with the
# WCA potential settings. In the last two lines
# we set up a WCA between the particles and the 
# walls

set sig 1.0
set cut [expr 1.12246*$sig]
set eps 1.0
set shift [expr 0.25*$eps]

inter 0 0 lennard-jones $eps $sig $cut $shift 0
inter 1 0 lennard-jones $eps $sig $cut $shift 0
inter 1 1 lennard-jones $eps $sig $cut $shift 0

inter 0 2 lennard-jones $eps $sig $cut $shift 0
inter 1 2 lennard-jones $eps $sig $cut $shift 0

# we require a electrostatic algorithm capable
# of handling the partial periodic boundary 
# conditions. In this tutorial we select the 
# MMM2D method

cellsystem layered 3
inter coulomb 1.0 mmm2d 1.0e-4

# we determine the number of degrees of freedom 
# available to the particle. This value is dependent
# on the incompiled features of ESPResSo

if { [regexp "ROTATION" [code_info]] } {
    set deg_free 6
} {
    set deg_free 3
}

# The warming of the system is far more complex for 
# the slit-pore system. We cannot simply ramp up the
# WCA interactions for all elements, since this would
# allow the particles to penetrate the walls and escape
# their confinement. We therefore set a capping value
# for the particles only, and a radius for this cap. 
# We also slowly increase the Bjerrum length for the 
# system. Finally, we kill the force cap and set the
# final value for the Bjerrum lengths.

set tcl_precision 8
inter forcecap individual
set integ_steps 200
puts "Warming loop"

for {set cap 1} {$cap < 10} {incr cap} {
  set rad [expr 1.0 - 0.5*$cap/10.0]
  set lb [expr 1.0*$cap/10.0]

  inter 0 0 lennard-jones $eps $sig $cut $shift 0 $rad
  inter 1 0 lennard-jones $eps $sig $cut $shift 0 $rad
  inter 1 1 lennard-jones $eps $sig $cut $shift 0 $rad
  inter coulomb $lb mmm2d 1e-4

  integrate $integ_steps

  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
  puts -nonewline "t = [setmd time], E = [analyze energy total], T = $temp \r"
  flush stdout
}

puts "\n"
inter forcecap 0
inter coulomb 1.0 mmm2d 1e-4

# next, we equilibrate the system that we
# obtained after the warming process

set integ_steps 200
puts "Equilibration loop"

for {set i 0} { $i < 16 } { incr i} {
  integrate $integ_steps
  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
  puts -nonewline "t = [setmd time], E = [analyze energy total], T = $temp \r"
  flush stdout
}

puts "\n"

# finally we determine physical quanties
# in a production run

set pro [open "data/walls.dat" "w"]
puts $pro "# Time / Total Energy / Kinetic Temperature"

set vtf_file [open "data/charged_walls.vtf" "w"]
writevsf $vtf_file ignore_charges

set integ_steps 20
puts "Production loop"

for {set i 0} { $i < 2000 } { incr i} {
  integrate $integ_steps

  set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
  puts $pro "[setmd time] [analyze energy total] $temp "

  puts -nonewline "t = [setmd time], E = [analyze energy total], T = $temp \r"
  flush stdout

  if {$i%10==0} {
    set f [open "data/config_2d_[expr $i/10]" "w"]
    blockfile $f write tclvariable {box_lxy box_lz density}
    blockfile $f write variable box_l
    blockfile $f write particles {id pos type}
    close $f

    writevcf $vtf_file folded
  }
}

close $pro
close $vtf_file
puts "\n"

exit
