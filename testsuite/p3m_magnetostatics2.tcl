# Copyright (C) 2010,2011 The ESPResSo project
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

require_feature "MAGNETOSTATICS" 
require_feature "FFTW"
require_feature "ROTATION"

set tcl_precision 15
  
# volume fraction
set rho 0.09
set dipole_lambda 3.0
set particle_radius 1
set bjerrum 1
set int_time_step 0.002

set accuracy_p3m 1.0e-4
set max_allowed_rel_error $accuracy_p3m

set pos_dip_data [open "p3m_magnetostatics2_system.data" r]
  
set counter 0
while { [eof $pos_dip_data] == 0 } {
    set line [gets $pos_dip_data]
    incr counter
}
close $pos_dip_data
set n_particle [expr $counter -1]

set dipole_modulus [expr sqrt($dipole_lambda * pow(2*$particle_radius,3))]

# setting box paramters
set box_l [expr pow(((4 * $n_particle * 3.141592654) / (3*$rho)), 1.0/3.0)*$particle_radius]  
set skin 0.5

# give Espresso some parameters
setmd time_step $int_time_step
setmd skin $skin
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1
setmd max_num_cells 2500

# read positions and torques from file
set pos_dip_data [open "p3m_magnetostatics2_system.data" r]
for {set i 0} {$i < $n_particle} {incr i} { 
    set line [gets $pos_dip_data] 
    set partid [lindex $line 0]
    set posx [lindex $line 1]  
    set posy [lindex $line 2] 
    set posz [lindex $line 3] 
    set posx [expr $posx + ($box_l / 2.0)] 
    set posy [expr $posy +($box_l /  2.0)]
    set posz [expr $posz +($box_l /  2.0)]
    
    set dipx  [lindex $line 4]
    set dipy  [lindex $line 5] 
    set dipz  [lindex $line 6] 
    part $partid pos $posx $posy $posz dip $dipx $dipy $dipz
} 
close $pos_dip_data 

thermostat off 

puts "Dipolar P3M test case is running ..."

#puts [inter magnetic $bjerrum p3m tunev2 mesh 32 accuracy $accuracy_p3m]
#puts [inter magnetic]

# the following parameters where obtained by running the above tuning
inter magnetic $bjerrum p3m 10.8941877883205 32 5 0.301164191913605 9.92326028926542e-05

inter magnetic epsilon metallic

# set kkk [inter magnetic] 
# set xxx [split $kkk {" "}]
# set C_r_cut [lindex $xxx 3]   
# set C_mesh  [lindex $xxx 4]  
# set C_cao   [lindex $xxx 5]  
# set C_alpha [lindex $xxx 6]

integrate 0

# initialize 
set avr_rel_err_fx_i 0
set avr_rel_err_fy_i 0
set avr_rel_err_fz_i 0

set avr_rel_err_tx_i 0
set avr_rel_err_ty_i 0
set avr_rel_err_tz_i 0

set avr_diff_fx_i 0
set avr_diff_fy_i 0
set avr_diff_fz_i 0

set avr_diff_tx_i 0
set avr_diff_ty_i 0
set avr_diff_tz_i 0

set failed 0


# differences of forces and torques wich are smaller
# than than min_diff are considered as zero
set min_diff [expr $accuracy_p3m * 10] 

# compare forces and torques with results from "Wang code" 
set force_torque_data [open "p3m_magnetostatics2_expected.data" r]
for {set i 0} {$i < $n_particle} {incr i} {
    set line [gets $force_torque_data]
    
    set rel_err_force_x [expr  ([lindex [part $i print force] 0] - [lindex $line 1])/ [lindex $line 1] ]
    set diff_force_x [expr sqrt(([lindex [part $i print force] 0] - [lindex $line 1])*([lindex [part $i print force] 0] - [lindex $line 1]))]
    
    set rel_err_force_y [expr  ([lindex [part $i print force] 1] - [lindex $line 2])/ [lindex $line 2] ] 
    set diff_force_y [expr sqrt(([lindex [part $i print force] 1] - [lindex $line 2])*([lindex [part $i print force] 1] - [lindex $line 2]))]   
    
    set rel_err_force_z [expr  ([lindex [part $i print force] 2] - [lindex $line 3])/ [lindex $line 3] ] 
    set diff_force_z [expr sqrt(([lindex [part $i print force] 2] - [lindex $line 3])*([lindex [part $i print force] 2] - [lindex $line 3]))]  
    
    
    if { ($rel_err_force_x > $max_allowed_rel_error & $diff_force_x > $min_diff ) || ($rel_err_force_y > $max_allowed_rel_error & $diff_force_y > $min_diff) || ($rel_err_force_z > $max_allowed_rel_error & $diff_force_z > $min_diff) } {
        puts "\n\nrelative error of FORCE acting on particle $i is bigger than $max_allowed_rel_error:"
	puts "         rel_err = ($rel_err_force_x, $rel_err_force_y, $rel_err_force_y )  "
        puts "         difference = ($diff_force_x, $diff_force_y, $diff_force_z)"
        set failed 1
    }
    
    set avr_rel_err_fx_i [expr $avr_rel_err_fx_i + sqrt($rel_err_force_x * $rel_err_force_x)]
    set avr_rel_err_fy_i [expr $avr_rel_err_fy_i + sqrt($rel_err_force_y * $rel_err_force_y)]
    set avr_rel_err_fz_i [expr $avr_rel_err_fz_i + sqrt($rel_err_force_z * $rel_err_force_z)]
    
    set avr_diff_fx_i [expr $avr_diff_fx_i + $diff_force_x] 
    set avr_diff_fy_i [expr $avr_diff_fy_i + $diff_force_y] 
    set avr_diff_fz_i [expr $avr_diff_fz_i + $diff_force_z] 
    
    set rel_err_torque_x [expr ([lindex [part $i print torque] 0] - [lindex $line 4])/[lindex $line 4] ]
    set diff_torque_x [expr sqrt(([lindex [part $i print torque] 0] - [lindex $line 4])*([lindex [part $i print torque] 0] - [lindex $line 4]))]

    set rel_err_torque_y [expr ([lindex [part $i print torque] 1] - [lindex $line 5])/[lindex $line 5] ]
    set diff_torque_y [expr sqrt(([lindex [part $i print torque] 1] - [lindex $line 5])*([lindex [part $i print torque] 1] - [lindex $line 5]))]
    
    set rel_err_torque_z [expr ([lindex [part $i print torque] 2] - [lindex $line 6])/[lindex $line 6] ]
    set diff_torque_z [expr sqrt(([lindex [part $i print torque] 2] - [lindex $line 6])*([lindex [part $i print torque] 2] - [lindex $line 6]))]
    
    if { ($rel_err_torque_x > $max_allowed_rel_error & $diff_torque_x > $min_diff) || ($rel_err_torque_y > $max_allowed_rel_error & $diff_torque_y >  $min_diff) || ($rel_err_torque_z > $max_allowed_rel_error& $diff_torque_z > $min_diff)} {
        puts stderr "\nrelative error of TORQUE acting on particle $i is bigger than $max_allowed_rel_error:"
	puts stderr "       rel_err = ($rel_err_torque_x, $rel_err_torque_y, $rel_err_torque_z)"
	puts stderr "       difference = ($diff_torque_x, $diff_torque_y, $diff_torque_z)"
        set failed 1 
    }
    
    set avr_rel_err_tx_i [expr $avr_rel_err_tx_i + sqrt($rel_err_torque_x * $rel_err_torque_x)]
    set avr_rel_err_ty_i [expr $avr_rel_err_ty_i + sqrt($rel_err_torque_y * $rel_err_torque_y)]
    set avr_rel_err_tz_i [expr $avr_rel_err_tz_i + sqrt($rel_err_torque_z * $rel_err_torque_z)]
    
    set avr_diff_tx_i [expr $avr_diff_tx_i + $diff_torque_x]
    set avr_diff_ty_i [expr $avr_diff_ty_i + $diff_torque_y]
    set avr_diff_tz_i [expr $avr_diff_tz_i + $diff_torque_z]
}
close $force_torque_data

if {$failed==1} {

    #set avr_rel_err_fx [expr ($avr_rel_err_fx_i / $n_particle )] 
    #set avr_rel_err_fy [expr ($avr_rel_err_fy_i / $n_particle )]
    #set avr_rel_err_fz [expr ($avr_rel_err_fz_i / $n_particle )]

    #puts stderr " "
    #puts stderr "avr rel_force_err x is:  $avr_rel_err_fx"
    #puts stderr "avr rel_force_err y is:  $avr_rel_err_fy"
    #puts stderr "avr rel_force_err z is:  $avr_rel_err_fz"

    #set avr_rel_err_tx [expr ($avr_rel_err_tx_i / $n_particle )]
    #set avr_rel_err_ty [expr ($avr_rel_err_ty_i / $n_particle )]
    #set avr_rel_err_tz [expr ($avr_rel_err_tz_i / $n_particle )]
    #puts stderr " "
    #puts stderr "avr rel_torque_err x is:  $avr_rel_err_tx"    
    #puts stderr "avr rel_torque_err y is:  $avr_rel_err_ty" 
    #puts stderr "avr rel_torque_err z is:  $avr_rel_err_tz" 
    #puts stderr " "

    #set avr_diff_fx [expr $avr_diff_fx_i/$n_particle]
    #set avr_diff_fy [expr $avr_diff_fy_i/$n_particle]
    #set avr_diff_fz [expr $avr_diff_fz_i/$n_particle]
    #puts stderr " "
    #puts stderr "avr difference force x is: $avr_diff_fx"
    #puts stderr "avr difference force y is: $avr_diff_fy"
    #puts stderr "avr difference force z is: $avr_diff_fz"

    #set avr_diff_tx [expr $avr_diff_tx_i/$n_particle]
    #set avr_diff_ty [expr $avr_diff_ty_i/$n_particle]
    #set avr_diff_tz [expr $avr_diff_tz_i/$n_particle]
    #puts " "
    #puts stderr "avr difference torque x is: $avr_diff_tx"
    #puts stderr "avr difference torque y is: $avr_diff_ty"
    #puts stderr "avr difference torque z is: $avr_diff_tz"

    puts stderr "\n\n\n\n****************************************************************************************"
    puts stderr "\n relative error of forces and torques was larger than $max_allowed_rel_error, the test FAILED !"
    puts stderr "\n see some details above "
    puts stderr "\n***************************************************************************************"

}

#if {$failed==0} { 
#    puts "\n\n*******************************************************************************************"
#    puts "\n  relative error of forces and torques was smaller than $max_allowed_rel_error, the test was SUCCESSEFUL"
#    puts "\n*******************************************************************************************"
#}
#
#puts "\n\nthis test was dones, for:"
#puts "\n     number of particles: $n_particle"
#puts "     accuracy of p3m was set to: $accuracy_p3m"
#     
#    puts "\ntuning parameters for magnetic P3M was set via autotuning (v2) to:" 
#    puts "\n       rcut=$C_r_cut"
#    puts "       mesh=$C_mesh"
#    puts "       cao=$C_cao"
#    puts "       alpha=$C_alpha"




#puts "\n\n\n"


