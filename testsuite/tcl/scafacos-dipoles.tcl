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
require_feature "ROTATION"
require_feature "SCAFACOS_DIPOLES"

set tcl_precision 15

# volume fraction
set rho 0.09
set dipole_lambda 3.0
set particle_radius 1.0
set bjerrum 1
set int_time_step 0.002

# this are NOT relative errors, note!
set accuracy_p3m 1.0e-4

set pos_dip_data [open "p3m_magnetostatics2_system.data" r]

set counter 0
while { [eof $pos_dip_data] == 0 } {
    set line [gets $pos_dip_data]
    incr counter
}
close $pos_dip_data
set n_particle [expr $counter-1]

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


inter magnetic $bjerrum scafacos p2nfft p2nfft_verbose_tuning 0 pnfft_N 32,32,32 pnfft_window_name bspline pnfft_m 5 p2nfft_ignore_tolerance 1 pnfft_diff_ik 1 p2nfft_r_cut 11 p2nfft_alpha 0.3

puts [time { integrate 0} ]

set failed 0

set rms_force  0
set rms_torque 0

if { [catch {

    # compare forces and torques with results from "Wang code" 
    set force_torque_data [open "p3m_magnetostatics2_expected.data" r]
    for {set i 0} {$i < $n_particle} {incr i} {
	set line [gets $force_torque_data]
	
	set diff_force [vecsub [part $i print force] [lrange $line 1 3]]
	set error_force [veclen $diff_force]
	set rms_force [expr $rms_force + pow($error_force, 2)]

	set diff_torque [vecsub [part $i print torque_lab] [lrange $line 4 6]]
	set error_torque [veclen $diff_torque]
	set rms_torque [expr $rms_torque + pow($error_torque, 2)]
    }
    close $force_torque_data

    set rms_force [expr sqrt($rms_force/$n_particle)]
    if {$rms_force > $accuracy_p3m} {
	puts "rms force error $rms_force -> too large"
	set failed 1
    }

    set rms_torque [expr sqrt($rms_torque/$n_particle)]
    if {$rms_torque > 2*$accuracy_p3m} {
	puts "rms torque error $rms_torque -> too large"
	set failed 1
    }

    set E [analyze energy magnetic]
    # Compare energy against result of p3m_magntostatics2.tcl
    set E_diff [expr abs($E-5.571)]
    if { $E_diff >1E-3 } {
      set failed 1
      puts "Energy difference too large" 
    }


    if {$failed==1} {
	error_exit "\nerror of forces, torques or magnetic energy was too large"
    }
} res ] } {
    error_exit $res
}

puts "rms force diff: $rms_force" 
puts "rms torque diff: $rms_torque" 
puts "energy diff: $E_diff" 

exit 0
