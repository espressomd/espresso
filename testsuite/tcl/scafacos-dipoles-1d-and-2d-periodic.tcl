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


# This tests the scafacos p2nfft dipolar calculations by matching against
# reference data from direct summation. In 2d, reference data from the mdlc
# test case is used




source "tests_common.tcl"

require_feature "DIPOLES" 
require_feature "FFTW"
require_feature SCAFACOS_DIPOLES
require_feature PARTIAL_PERIODIC
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
setmd max_num_cells 2500


foreach dim {1 2} {

puts "${dim}d periodic"

# Read reference data
if { $dim == 2} {
  set file_prefix "mdlc"
  setmd periodic 1 1 0
} else {
  setmd periodic 1 0 0
  set file_prefix "scafacos_dipoles_1d"
}
set f [open ${file_prefix}_reference_data_energy.dat]
set ref_E [gets $f]
close $f


# Particles
set ids ""
set f [open "${file_prefix}_reference_data_forces_torques.dat"]
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
#inter magnetic 1 p3m 1.5 128 6 2.855458 8.437322E-6
if {$dim ==  2} {
 inter magnetic 1 scafacos p2nfft p2nfft_verbose_tuning 0 pnfft_N 64,64,160 pnfft_window_name bspline pnfft_m 6 p2nfft_ignore_tolerance 1 pnfft_diff_ik 0 p2nfft_r_cut 6 p2nfft_alpha 1.5 p2nfft_epsB 0.1
} else {
  if { $dim == 1} {
    # 1d periodic in x
    inter magnetic 1 scafacos p2nfft p2nfft_verbose_tuning 1 pnfft_N 32,128,128 pnfft_direct 0 p2nfft_r_cut 2.855 p2nfft_alpha 1.5 p2nfft_intpol_order -1 p2nfft_reg_kernel_name ewald p2nfft_p 16 p2nfft_ignore_tolerance 1 pnfft_window_name bspline pnfft_m 8 pnfft_diff_ik 1 p2nfft_epsB 0.125
 } else {
   puts "Wrong dimensions $dim"
   exit 1
 }
}
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
inter magnetic 0 
part delete
}





