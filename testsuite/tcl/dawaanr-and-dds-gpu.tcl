

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

proc stopAll {} {
  for {set i 0} {$i <=[setmd max_part] } {incr i} {
    if { [part $i print pos] != "na" } {
      part $i v 0 0 0 omega 0 0 0
    }
  }
}

proc error { msg } {
  puts $msg
}

source "tests_common.tcl"
require_feature "DIPOLES"
require_feature "CUDA"
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
puts "- Testcase dawaanr-and-dds-gpu.tcl"
puts "----------------------------------------------"

set pf_dds_gpu 2.34
set pf_dawaanr 3.524
set ratio_dawaanr_dds_gpu [expr $pf_dawaanr / $pf_dds_gpu]

set l 15 
setmd box_l $l $l $l
setmd periodic 0 0 0
setmd time_step 0.0001
setmd skin 0.1


foreach n { 110 111 540 541 } {



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
minimize_energy 0 500 0.1 0.1
stopAll

inter 0 0 lennard-jones 0 0 0 0



setmd skin 0
setmd time_step 0.01

inter magnetic $pf_dawaanr dawaanr
integrate 0 recalc_forces

set dawaanr_f ""
set dawaanr_t ""

for {set i 0} {$i<$n} {incr i} {
  lappend dawaanr_f [part [expr 2*$i] print f]
  lappend dawaanr_t [part [expr 2*$i] print torque]
}
set dawaanr_e [analyze energy total]

inter magnetic 0

integrate 0 recalc_forces

inter magnetic $pf_dds_gpu dds-gpu

integrate 0 recalc_forces

set ddsgpu_f ""
set ddsgpu_t ""

for {set i 0} {$i<$n} {incr i} {
  lappend ddsgpu_f [part [expr 2*$i] print f]
  lappend ddsgpu_t [part [expr 2*$i] print torque]
}

set ddsgpu_e [analyze energy total]

# compare

for {set i 0} {$i<$n} {incr i} {
  if { ! [vectorsTheSame [lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_dds_gpu [lindex $ddsgpu_f $i]]] } {
    error "Forces on particle don't match. $i  [vecsub [lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_dds_gpu [lindex $ddsgpu_f $i]]]\n[lindex $dawaanr_f $i] [vecscale $ratio_dawaanr_dds_gpu [lindex $ddsgpu_f $i]]"
  }
  if { ! [vectorsTheSame [lindex $dawaanr_t $i] [vecscale $ratio_dawaanr_dds_gpu [lindex $ddsgpu_t $i]]] } {
    error "torques on particle $i don't match.  [vecsub [lindex $dawaanr_t $i] [lindex $ddsgpu_t $i]] "
  }
}

if { abs($dawaanr_e - $ddsgpu_e*$ratio_dawaanr_dds_gpu) > 0.001 } {
  error "Energies for dawaanr $dawaanr_e and dds_gpu $ddsgpu_e don't match."
}

# Does the energy stay constant when swapping particles
for {set i 0} {$i<1000} {incr i} {
  set a [expr int(rand() *$n)*2]
  set b [expr int(rand() *$n)*2]
  set dipa [part $a print dip]
  set posa [part $a print pos]
  eval "part $a pos [part $b print pos] dip [part $b print dip]"
  eval "part $b pos $posa dip $dipa"
  integrate 0 recalc_forces
  set E [analyze energy total]
  if { abs($E -$ddsgpu_e) > 1E-3 } {
    error "energy mismatch: $E $ddsgpu_e"
  }
}


integrate 0 recalc_forces
inter magnetic 0 
part deleteall


}
