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

source "tests_common.tcl"

require_feature "LEES_EDWARDS"
require_feature "LENNARD_JONES"
require_feature "DPD"

puts "------------------------------------------"
puts "- Testcase lees_edwards_dpd.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------"


# system setup
set L                 30     ;#box size
set temperature        1.0   ;#temperature in reduced units
set T                500     ;#temperature in Kelvin
set n_part          1400     ;#number of particles

set k_temp [expr $T * 0.00198721]
set time_step              0.001
set skin                   0.2

#############################################################
# Energy Parameters                                         #
#############################################################

# Epsilon in KCal/mol divided by KbT in KCal/mol
set lj_epsilon [expr  0.112/$k_temp]    

#############################################################
#  Lennard Jones Parameters                                 #
#############################################################
set lj1_sig 1.0
set lj1_cut [expr 2.5*$lj1_sig]
set lj1_shift "auto"

# set up global parameters
setmd box_l     $L $L $L
setmd time_step $time_step
setmd skin      $skin

# Seed the Random Number Generator
set ran_seed 1902
set cmd "t_random seed"
for {set i 0} {$i < [setmd n_nodes]} { incr i } { 
   lappend cmd [expr $ran_seed + $i] 
}
eval $cmd


#############################################################
#  Interaction Setup                                        #
#############################################################

inter 0 0 lennard-jones $lj_epsilon $lj1_sig $lj1_cut auto 0  

# Place particles randomly
for { set i 0 } { $i < $n_part } { incr i } {
    part $i pos [expr [t_random]*$L] [expr [t_random]*$L] [expr [t_random]*$L]
    part $i type 0
}

##################################################
# simulation

set max_step_shear     20000
set shear_equil         5000
set mean_from           5000
set shear_per              1
set write_per            100
set shear_rate             0.005
set offset                 0.0
set cap_until           5000

#modified Langevin thermostat assumes a linear flow profile
thermostat langevin $temperature 5.0
puts "#Shear displacement   |  <vx^2>    <vy^2>   <vz^2>"
set fCap 0
for { set step 0 } { $step < $shear_equil } { incr step $shear_per } {
    if { $step % 200 == 0 && $fCap < 500 } then { 
       incr fCap 10
       inter forcecap  $fCap
    }

    lees_edwards_offset  $offset
    integrate            $shear_per

    if { [expr $step % $write_per] == 0 } then {
        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set vel [ part $i print v ] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]

        puts "$offset     $KEx $KEy $KEz"
        
    }

    set offset [expr $offset + $shear_rate * $shear_per]
}

##Below is the loop for an extended test, use this to
##see if the code is stable without a forcecap for a long time.

##
thermostat off
thermostat dpd $temperature 10.0 2.5

puts "#Shear displacement   |  <vx^2>    <vy^2>   <vz^2>"
set mean_x2 0.0
set mean_y2 0.0
set mean_z2 0.0
set count   0.0

for { set step 0 } { $step < $max_step_shear } { incr step $shear_per } {
    lees_edwards_offset  $offset
    integrate            $shear_per

    if { $step == $cap_until } then { 
	puts "Removing forcecap"
        inter forcecap 0
    }

    if { [expr $step % $write_per] == 0 } then {
        set vx2 0.0
        set vy2 0.0
        set vz2 0.0
        for { set i 0 } { $i < $n_part } { incr i } {
            set vel [ part $i print v ] 
            set vx2 [expr $vx2 + [expr [lindex $vel 0] * [lindex $vel 0]]] 
            set vy2 [expr $vy2 + [expr [lindex $vel 1] * [lindex $vel 1]]] 
            set vz2 [expr $vz2 + [expr [lindex $vel 2] * [lindex $vel 2]]] 
        }

        set KEx             [expr 2.0 * $vx2 / (3 * $n_part )]
        set KEy             [expr 2.0 * $vy2 / (3 * $n_part )]
        set KEz             [expr 2.0 * $vz2 / (3 * $n_part )]

        if {$step > $cap_until} {
            set mean_x2 [expr $mean_x2 + $KEx ]
            set mean_y2 [expr $mean_y2 + $KEy ]
            set mean_z2 [expr $mean_z2 + $KEz ]
            set count   [expr $count + 1.0 ]

            set mmX [expr $mean_x2 / $count]
        } {
            set mmX "not yet computed"
        }
        puts "$offset     $KEx $KEy $KEz running mean vx^2:  $mmX"
    }
    set offset [expr $offset + $shear_rate * $shear_per]
}

set mean_x2 [expr $mean_x2 / $count ] 
set mean_y2 [expr $mean_y2 / $count ] 
set mean_z2 [expr $mean_z2 / $count ] 
puts "#means:   $mean_x2 $mean_y2 $mean_z2"
puts "#refvals: 2.15(pm 0.1) 0.68(pm 0.05) 0.68(pm 0.05)"

if { abs( $mean_x2 - 2.15) > 0.1 } {
    error_exit "Lees Edwards Test failed"
}
if { abs( $mean_y2 - 0.68) > 0.05 } {
    error_exit "Lees Edwards Test failed"
}
if { abs( $mean_z2 - 0.68) > 0.05 } {
    error_exit "Lees Edwards Test failed"
}

exit 0
