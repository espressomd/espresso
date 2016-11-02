# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

source "tests_common.tcl"

require_feature "DIPOLES" 
require_feature "FFTW"
require_feature "ROTATION"
require_feature "CONSTRAINTS"
if {[has_feature "LEES_EDWARDS"]} {
    require_max_nodes_per_side 1
} {
    require_max_nodes_per_side 2
}

proc vectorsTheSame {a b} {
 set tol 3E-3
 set diff [vecsub $a $b]
 if { [veclen $diff] > $tol } {
  return 0
  puts "Difference: [veclen $diff]"
 }
 return 1
}
set tcl_precision 15



##############################################
# volume fraction
set rho 0.3

# This is only for box size calculation. The actual particle numbwe is
# lower, because particles are removed from the mdlc gap region
set n_particle 100

set particle_radius 0.5
set dipole_lambda 3.0

#################################################

set box_l [expr pow(((4 * $n_particle * 3.141592654) / (3*$rho)), 1.0/3.0)*$particle_radius]  
puts "box_lenght $box_l"
set skin 0.5


# give Espresso some parameters
setmd time_step 0.01
setmd skin $skin
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1
setmd max_num_cells 2500



# Read reference data
set f [open mdlc_reference_data_energy.dat]
set ref_E [gets $f]
close $f


# Particles
set ids ""
set f [open "mdlc_reference_data_forces_torques.dat"]
while { ! [eof $f] } {
  set line [gets $f]
  set id [lindex $line 0]
  if {$id ==""} {
    continue
  }
  set pos [lrange $line 1 3]
  set dip [lrange $line 4 6]
  set forces($id) [lrange $line 7 9]
  set torques($id) [lrange $line 10 12]
  set cmd "part $id pos $pos dip $dip"
  eval $cmd
  lappend ids $id

}
#puts [inter magnetic 1 p3m tunev2 accuracy 1E-5 mesh 128 r_cut 1.5 ]
inter magnetic 1 p3m 1.5 64 7 2.855458 8.437322E-6
inter magnetic epsilon metallic
puts [inter magnetic mdlc 1E-7 2]
puts [inter magnetic]
#inter magnetic 1 mdds n_cut 20 
thermostat off
integrate 0

set err_f 0.
set err_t 0.
foreach i $ids {
 set err_f [expr $err_f +[veclen [vecsub [part $i print force] $forces($i)]]]
 set err_t [expr $err_f +[veclen [vecsub [part $i print torque_lab] $torques($i)]]]
}
set err_f [expr $err_f /sqrt([setmd n_part])]
set err_t [expr $err_f /sqrt([setmd n_part])]
set err_e [expr abs([analyze energy magnetic] - $ref_E)]
puts "Force error: $err_f"
puts "torque error: $err_t"
puts "energy error: $err_e"

# Tolerance values
set tol_f 1E-3
set tol_t 1E-3
set tol_e 1E-3

if { $err_f > $tol_f } {
 error "Force error too large"
}

if { $err_t > $tol_t } {
 error "Torque error too large"
}
if { $err_e > $tol_e } {
 error "Energy error too large"
}






