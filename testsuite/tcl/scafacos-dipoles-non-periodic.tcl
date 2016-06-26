

# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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


# This tests scafacos p2nfft dipolar calculations with open boundaries
# by matching against the direct summation (dawaanr)



proc stopAll {} {
  for {set i 0} {$i <=[setmd max_part] } {incr i} {
    if { [part $i print pos] != "na" } {
      part $i v 0 0 0 omega 0 0 0
    }
  }
}

proc error { msg } {
  error_exit $msg
}

source "tests_common.tcl"
require_feature "DIPOLES"
require_feature "SCAFACOS_DIPOLES"
require_feature "PARTIAL_PERIODIC"
require_max_nodes_per_side 1

source "tests_common.tcl"

set tcl_precision 14

proc vectorsTheSame {a b} {
 set tol 3E-3
 set diff [vecsub $a $b]
 if { [veclen $diff] > $tol } {
  return 0
  puts "Difference: [veclen $diff]"
 }
 return 1
}


puts "----------------------------------------------"
puts "- Testcase scafacos-dipoles-non-periodic "
puts "----------------------------------------------"

#set pf_scafacos 2.34
#set pf_dawaanr 3.524
set pf_scafacos 1.
set pf_dawaanr 1.
set ratio_dawaanr_scafacos [expr $pf_dawaanr / $pf_scafacos]

set l 15 
setmd box_l $l $l $l
setmd periodic 0 0 0
setmd time_step 0.0001
setmd skin 0.1


foreach n { 110 541 } {



puts "$n particles "
set dipole_modulus 1.3
for {set i 0 } {$i < $n} {incr i} {
 set posx [expr $l*[expr [t_random]]]
 set posy [expr $l*[expr [t_random]]]
 set posz [expr $l*[expr [t_random]]]
 set costheta [expr 2*[expr [t_random]] - 1]
 set sintheta [expr sin(acos($costheta))]
 set phi [expr 2*3.1415926536*[expr [t_random]]] 
 set dipx [expr $sintheta*cos($phi)*$dipole_modulus]
 set dipy [expr $sintheta*sin($phi)*$dipole_modulus]
 set dipz [expr $costheta*$dipole_modulus]

 part [expr 2*$i] pos $posx $posy $posz type 0 dip $dipx $dipy $dipz v 0 0 0
}
inter 0 0 lennard-jones 10 0.5 0.55 auto
thermostat langevin 0 10 
setmd periodic 1 1 1
minimize_energy 0 500 0.1 0.1

stopAll
# Fold particles back into box
for {set i 0 } {$i < $n} {incr i} {
 eval "part [expr 2*$i] pos [part [expr 2*$i] print folded_pos]" 
}
setmd periodic 0 0 0

inter 0 0 lennard-jones 0 0 0 0



setmd skin 0
setmd time_step 0.01

inter magnetic $pf_dawaanr dawaanr
puts "Dawaanr [time {integrate 0 recalc_forces}]
"
set dawaanr_f ""
set dawaanr_t ""

for {set i 0} {$i<$n} {incr i} {
  lappend dawaanr_f [part [expr 2*$i] print f]
  lappend dawaanr_t [part [expr 2*$i] print torque]
}
set dawaanr_e [analyze energy total]

inter magnetic 0

integrate 0 recalc_forces
# Obtain reference parameters from p3m
#setmd periodic 1 1 1
#puts [inter magnetic 1 p3m tunev2 accuracy 1E-5 mesh 64 ]
#exit
inter magnetic $pf_scafacos scafacos p2nfft p2nfft_verbose_tuning 0 pnfft_N 128,128,128 pnfft_window_name bspline pnfft_m 7 p2nfft_ignore_tolerance 1 p2nfft_r_cut 5 p2nfft_alpha 1 p2nfft_p 8 p2nfft_epsB 0.0625 

puts "p2nfft [time {integrate 0 recalc_forces}]"

set scafacos_f ""
set scafacos_t ""

for {set i 0} {$i<$n} {incr i} {
  lappend scafacos_f [part [expr 2*$i] print f]
  lappend scafacos_t [part [expr 2*$i] print torque]
}

set scafacos_e [analyze energy total]

# compare

for {set i 0} {$i<$n} {incr i} {
#   puts "[expr [veclen [lindex $dawaanr_t $i]] / ($ratio_dawaanr_scafacos * [veclen [lindex $scafacos_t $i]])]"
#  puts "[lindex $dawaanr_f $i] | [vecscale $ratio_dawaanr_scafacos [lindex $scafacos_f $i]]"
  if { ! [vectorsTheSame [lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_scafacos [lindex $scafacos_f $i]]] } {
    error "Forces on particle don't match. $i  [vecsub [lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_scafacos [lindex $scafacos_f $i]]]\n[lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_scafacos [lindex $scafacos_f $i]]"
  }
  if { ! [vectorsTheSame [lindex $dawaanr_t $i] [vecscale $ratio_dawaanr_scafacos [lindex $scafacos_t $i]]] } {
    error "torques on particle $i don't match.  [vecsub [lindex $dawaanr_t $i] [lindex $scafacos_t $i]] "
  }
}

if { abs($dawaanr_e - $scafacos_e*$ratio_dawaanr_scafacos) > 0.001 } {
  error "Energies for dawaanr $dawaanr_e and scafacos $scafacos_e don't match."
}


# Does the energy stay constant when swapping particles
for {set i 0} {$i<10} {incr i} {
  set a [expr int([t_random] *$n)*2]
  set b [expr int([t_random] *$n)*2]
  set dipa [part $a print dip]
  set posa [part $a print pos]
  eval "part $a pos [part $b print pos] dip [part $b print dip]"
  eval "part $b pos $posa dip $dipa"
  integrate 0 recalc_forces
  set E [analyze energy total]
  if { abs($E -$scafacos_e) > 1E-3 } {
    error "energy mismatch: $E $scafacos_e"
  }
}


inter magnetic 0 
part deleteall


}
