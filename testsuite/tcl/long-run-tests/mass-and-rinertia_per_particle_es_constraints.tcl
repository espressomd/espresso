# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
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

require_feature "MASS"
require_feature "ROTATIONAL_INERTIA"
require_feature "LANGEVIN_PER_PARTICLE"

puts "------------------------------------------------"
puts "- Testcase mass-and-rinertia_per_particle_es_constraints.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"
puts "------------------------------------------------"

proc test_mass-and-rinertia_per_particle {test_case} {
    set PI 3.14159
    # make real random draw
    set cmd "t_random seed"
    for {set i 0} {$i < [setmd n_nodes]} { incr i } {
	  lappend cmd [expr [pid] + $i] }
    eval $cmd

    set gamma0 1.
    set gamma1 1. 
    
    # Decelleration
    setmd skin 0
    setmd time_step 0.007
    thermostat langevin 0 $gamma0 
    cellsystem nsquare 
    #-no_verlet_lists
    set J "10 10 10"
    set mass 12.74

    part deleteall
    part 0 pos 0 0 0 mass $mass rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1 v 1 1 1
    part 1 pos 0 0 0 mass $mass rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1 v 1 1 1
    puts "\n"
    switch $test_case {
        0 {
            puts "------------------------------------------------"
            puts "Test $test_case: no particle specific values"
            puts "------------------------------------------------"
            set gamma1 $gamma0 
        }
        1 {
            puts "------------------------------------------------"
            puts "Test $test_case: particle specific gamma but not temperature"
            puts "------------------------------------------------"
            part 0 gamma $gamma0
            part 1 gamma $gamma1
        }
        2 {
            puts "------------------------------------------------"
            puts "Test $test_case: particle specific temperature but not gamma"
            puts "------------------------------------------------"
            part 0 temp 0.0
            part 1 temp 0.0
            set gamma1 $gamma0 
        }
        3 {
            puts "------------------------------------------------"
            puts "Test $test_case: both particle specific gamma and temperature"
            puts "------------------------------------------------"
            part 0 gamma $gamma0 temp 0.0
            part 1 gamma $gamma1 temp 0.0
        }
    }
    setmd time 0
    for {set i 0} {$i <100} {incr i} {
        for {set k 0} {$k <3} {incr k} {
            if {[has_feature "VERLET_STEP4_VELOCITY"]} {
                set tolerance_o 1.25E-4
            } else {
                set tolerance_o 1.0E-2
            }
            set do0 [expr abs([lindex [part 0 print omega_body] $k] -exp(-$gamma0*[setmd time] /[lindex $J $k]))]
            set do1 [expr abs([lindex [part 1 print omega_body] $k] -exp(-$gamma1*[setmd time] /[lindex $J $k]))]
            if { $do0 > $tolerance_o || $do1 > $tolerance_o } {
                error_exit "Friction Deviation in omega too large. $i $k $do0 $do1"
            }
            if {[has_feature "VERLET_STEP4_VELOCITY"]} {
                set tolerance_v 4.5E-5
            } else {
                set tolerance_v 1.0E-2
            }
            set dv0 [expr abs([lindex [part 0 print v] $k] -exp(-$gamma0*[setmd time] / $mass))]
            set dv1 [expr abs([lindex [part 1 print v] $k] -exp(-$gamma1*[setmd time] / $mass))]
            if { $dv0 > $tolerance_v || $dv1 > $tolerance_v } {
                error_exit "Friction Deviation in v too large. $i $k $dv0 $dv1"
            }
        }
        integrate 10
    }


    #Accelerated motion
    # Left out, since not use of thermostat
    part deleteall


    # thermalization
    # Checks if every degree of freedom has 1/2 kT of energy, even when
    # mass and inertia tensor are active

    # 2 different langevin parameters for particles
    set gamma(0) [expr (0.2 + [t_random]) * 20]
    set gamma(1) [expr (0.2 + [t_random]) * 20]
    set temp(0) 2.5
    set temp(1) 2.0
    set gamma_rot_1(0) [expr (0.2 + [t_random]) * 20]
    set gamma_rot_2(0) [expr (0.2 + [t_random]) * 20]
    set gamma_rot_3(0) [expr (0.2 + [t_random]) * 20]
    set gamma_rot_1(1) [expr (0.2 + [t_random]) * 20]
    set gamma_rot_2(1) [expr (0.2 + [t_random]) * 20]
    set gamma_rot_3(1) [expr (0.2 + [t_random]) * 20]

    set box 35
    setmd box_l $box $box $box
    set kT 1.5
    if {[has_feature "ELECTROSTATICS"] && [has_feature "CONSTRAINTS"]} {
        cellsystem layered [expr 25/[setmd n_nodes]]
        setmd periodic 1 1 0
    }

    for {set k 0} {$k<2} {incr k} {
        set halfkT($k) 0.75
        if {$test_case == 2 || $test_case == 3} {
            set halfkT($k) [expr $temp($k)/2.]
        }        
        
        if {$test_case == 1 || $test_case == 3} {
            set gamma_tr($k) $gamma($k)
        } else {
            set gamma_tr($k) 1
        }
        # translational diffusion
        set D_tr($k) [expr 2 * $halfkT($k) / $gamma_tr($k)]
    }
    
    thermostat langevin $kT 1

    # no need to rebuild Verlet lists, avoid it
    setmd skin 1.0
    setmd time_step 0.004

    set n 200
    set mass [expr (0.2 + [t_random]) * 7]
    set j1 [expr (0.2 + [t_random]) * 7]
    set j2 [expr (0.2 + [t_random]) * 7]
    set j3 [expr (0.2 + [t_random]) * 7]

    for {set i 0} {$i<$n} {incr i} {
        for {set k 0} {$k<2} {incr k} {
            set ind [expr $i + $k*$n]
            part [expr $i + $k*$n] pos [expr [t_random] *$box] [expr [t_random] * $box] [expr $box * (0.15 + 0.7 * [t_random])] rinertia $j1 $j2 $j3 mass $mass omega_body 0 0 0 v 0 0 0
            switch $test_case {
                1 {part [expr $i + $k*$n] gamma $gamma($k) gamma_rot $gamma_rot_1($k) $gamma_rot_2($k) $gamma_rot_3($k)}
                2 {part [expr $i + $k*$n] temp $temp($k)}
                3 {part [expr $i + $k*$n] gamma $gamma($k) gamma_rot $gamma_rot_1($k) $gamma_rot_2($k) $gamma_rot_3($k) temp $temp($k)}
            }
            if {[has_feature "ELECTROSTATICS"]} {
                # Charges are assigned only to non-zero type of particles
                # part [expr $i + $k*$n] type $k q $k
                part [expr $i + $k*$n] type $k q 1
            }
        }
    }
    
    if {[has_feature "ELECTROSTATICS"] && [has_feature "CONSTRAINTS"]} {
        inter coulomb 1.0 mmm2d 1E-1
        constraint plate height [expr 1.1*$box] sigma [expr 1/(2.0*$PI)]
        
        set sig [expr 1.0];
        constraint wall normal 0 0 1 dist [expr 0.03 * $box] type 2 penetrable 0 reflecting 1
        constraint wall normal 0 0 -1 dist [expr -0.97 * $box] type 2 penetrable 0 reflecting 1
        
        set sig [expr 5.0/6.0]; set cut [expr 0.4*1.12246*$sig]
        set eps 5; set shift [expr 0.25*$eps]
        inter 2 1 lennard-jones $eps $sig $cut $shift 0
        inter 2 0 lennard-jones $eps $sig $cut $shift 0
    }
    
    prepare_vmd_connection "vmdfile" 10000

    for {set k 0} {$k<2} {incr k} {
        set vx2($k) 0.
        set vy2($k) 0.
        set vz2($k) 0.
        set ox2($k) 0.
        set oy2($k) 0.
        set oz2($k) 0.
        set dx2($k) 0.
        set dy2($k) 0.
        set dz2($k) 0.
        set dr($k) 0.
    }

    for {set p 0} {$p <$n} {incr p} {
        for {set k 0} {$k<2} {incr k} {
            set ind [expr $p + $k*$n]
            set pos0($ind) [part $ind print pos]
        }
    }
    set loops 1E3

    puts "Thermalizing..."
    set therm_steps 1200
    integrate $therm_steps
    puts "Measuring..."

    set int_steps 100
    for {set i 0} {$i <$loops} {incr i} {
        integrate $int_steps
        imd positions
        # Get kinetic energy in each degree of freedom for all particles
        for {set p 0} {$p <$n} {incr p} {
            for {set k 0} {$k<2} {incr k} {
                set ind [expr $p + $k*$n]
                set v [part [expr $p + $k*$n] print v]
                set o [part [expr $p + $k*$n] print omega_body]
                set pos [part [expr $p + $k*$n] print pos]
                set ox2($k) [expr $ox2($k) +pow([lindex $o 0],2)]
                set oy2($k) [expr $oy2($k) +pow([lindex $o 1],2)]
                set oz2($k) [expr $oz2($k) +pow([lindex $o 2],2)]
                set vx2($k) [expr $vx2($k) +pow([lindex $v 0],2)]
                set vy2($k) [expr $vy2($k) +pow([lindex $v 1],2)]
                set vz2($k) [expr $vz2($k) +pow([lindex $v 2],2)]
                set dx2($k) [expr pow([expr [lindex $pos 0] - [lindex $pos0($ind) 0]], 2)]
                set dy2($k) [expr pow([expr [lindex $pos 1] - [lindex $pos0($ind) 1]], 2)]
                set dz2($k) [expr pow([expr [lindex $pos 2] - [lindex $pos0($ind) 2]], 2)]
                set dt0 [expr $mass / $gamma_tr($k)]
                set dt [expr ($int_steps * ($i + 1) + $therm_steps) * [setmd time_step]]
                # translational diffusion variance: after a closed-form integration of the Langevin EOM
                if {!([has_feature "ELECTROSTATICS"] && [has_feature "CONSTRAINTS"])} {
                    set sigma2_tr($k) [expr $D_tr($k) * (6 * $dt + 3 * $dt0 * (-3 + 4 * exp(-$dt / $dt0) - exp(-2 * $dt / $dt0)))]
                    set dr($k) [expr $dr($k) + (($dx2($k)+$dy2($k)+$dz2($k)) - $sigma2_tr($k)) / $sigma2_tr($k)]
                } else {
                    if { $k == 0} {
                        set sigma2_tr($k) [expr $D_tr($k) * (2 * $dt + $dt0 * (-3 + 4 * exp(-$dt / $dt0) - exp(-2 * $dt / $dt0)))]
                        #set dr($k) [expr $dr($k) + ($dz2($k) - $sigma2_tr($k)) / $sigma2_tr($k)]
                        set dr($k) 0
                    } else {
                        # nothing to test due to contraints and the periodicity
                        set dr($k) 0
                    }
                }
            }
        }
    }

    set tolerance 0.15
    for {set k 0} {$k<2} {incr k} {
        set Evx($k) [expr 0.5 * $mass *$vx2($k)/$n/$loops]
        set Evy($k) [expr 0.5 * $mass *$vy2($k)/$n/$loops]
        set Evz($k) [expr 0.5 * $mass *$vz2($k)/$n/$loops]

        set Eox($k) [expr 0.5 * $j1 *$ox2($k)/$n/$loops]
        set Eoy($k) [expr 0.5 * $j2 *$oy2($k)/$n/$loops]
        set Eoz($k) [expr 0.5 * $j3 *$oz2($k)/$n/$loops]

        set dv($k) [expr 1./3. *($Evx($k) +$Evy($k) +$Evz($k))/$halfkT($k)-1.]
        set do($k) [expr 1./3. *($Eox($k) +$Eoy($k) +$Eoz($k))/$halfkT($k)-1.]
        set dr($k) [expr $dr($k)/$n/$loops]
        
        set dox($k) [expr ($Eox($k))/$halfkT($k)-1.]
        set doy($k) [expr ($Eoy($k))/$halfkT($k)-1.]
        set doz($k) [expr ($Eoz($k))/$halfkT($k)-1.]
        
        puts "\n"
        puts "1/2 kT = $halfkT($k)"
        puts "translation: $Evx($k) $Evy($k) $Evz($k) rotation: $Eox($k) $Eoy($k) $Eoz($k)"

        puts "Deviation in translational energy: $dv($k)"
        puts "Deviation in rotational energy: $do($k)"
        puts "Deviation in rotational energy per degrees of freedom: $dox($k) $doy($k) $doz($k)"
        puts "Deviation in translational diffusion: $dr($k) $mass $gamma_tr($k)"

        if { abs($dv($k)) > $tolerance } {
           error "Relative deviation in translational energy too large: $dv($k)"
        }
        if { abs($do($k)) > $tolerance } {
           puts "Moment of inertia principal components: $j1 $j2 $j3"
           error "Relative deviation in rotational energy too large: $do($k)"
        }
        if { abs($dr($k)) > $tolerance } {
           error "Relative deviation in translational dispersion too large: $dr($k)"
       }
       
        # SEMI_INTEGRATED is consistent for isotropic particles only
        if {![has_feature "SEMI_INTEGRATED"]} {
            if { abs($dox($k)) > $tolerance } {
            puts "Moment of inertia principal components: $j1 $j2 $j3"
            error "Relative deviation in rotational energy per the body axis X is too large: $dox($k)"
            }
            if { abs($doy($k)) > $tolerance } {
            puts "Moment of inertia principal components: $j1 $j2 $j3"
            error "Relative deviation in rotational energy per the body axis Y is too large: $doy($k)"
            }
            if { abs($doz($k)) > $tolerance } {
            puts "Moment of inertia principal components: $j1 $j2 $j3"
            error "Relative deviation in rotational energy per the body axis Z is too large: $doz($k)"
            }
        }
    }
}

# the actual testing
for {set i 0} {$i < 4} {incr i} {
    test_mass-and-rinertia_per_particle $i
}

exit 0
