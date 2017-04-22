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
require_feature "LANGEVIN_PER_PARTICLE"

puts "------------------------------------------------"
puts "- Testcase mass-and-rinertia_per_particle.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"
puts "------------------------------------------------"

proc test_mass-and-rinertia_per_particle {test_case} {
    # make real random draw
    set cmd "t_random seed"
    for {set i 0} {$i < [setmd n_nodes]} { incr i } {
	  lappend cmd [expr [pid] + $i] }
    eval $cmd

    set gamma(0) 1.
    set gamma(1) 1. 
    
    # Decelleration
    setmd skin 0
    setmd time_step 0.007
    thermostat langevin 0 $gamma(0)
    cellsystem nsquare 
    #-no_verlet_lists
    set J "10 10 10"
    set mass 12.74

    part deleteall
    if {[has_feature "ROTATIONAL_INERTIA"]} {
        part 0 pos 0 0 0 mass $mass rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1 v 1 1 1
        part 1 pos 0 0 0 mass $mass rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1 v 1 1 1
    } else {
        if {[has_feature "ROTATION"]} {
            part 0 pos 0 0 0 mass $mass rinertia [lindex $J 0] omega_body 1 1 1 v 1 1 1
            part 1 pos 0 0 0 mass $mass rinertia [lindex $J 0] omega_body 1 1 1 v 1 1 1
        } else {
            part 0 pos 0 0 0 mass $mass rinertia [lindex $J 0] v 1 1 1
            part 1 pos 0 0 0 mass $mass rinertia [lindex $J 0] v 1 1 1
        }
    }
    puts "\n"
    switch $test_case {
        0 {
            puts "------------------------------------------------"
            puts "Test $test_case: no particle specific values"
            puts "------------------------------------------------"
            set gamma1 $gamma(0)
        }
        1 {
            puts "------------------------------------------------"
            puts "Test $test_case: particle specific gamma but not temperature"
            puts "------------------------------------------------"
            part 0 gamma $gamma(0)
            part 1 gamma $gamma(1)
        }
        2 {
            puts "------------------------------------------------"
            puts "Test $test_case: particle specific temperature but not gamma"
            puts "------------------------------------------------"
            part 0 temp 0.0
            part 1 temp 0.0
            set gamma1 $gamma(0)
        }
        3 {
            puts "------------------------------------------------"
            puts "Test $test_case: both particle specific gamma and temperature"
            puts "------------------------------------------------"
            if {[has_feature "PARTICLE_ANISOTROPY"]} {
                part 0 gamma $gamma(0) $gamma(0) $gamma(0) temp 0.0
                part 1 gamma $gamma(1) $gamma(1) $gamma(1) temp 0.0
            } else {
                part 0 gamma $gamma(0) temp 0.0
                part 1 gamma $gamma(1) temp 0.0
            }
        }
    }
    setmd time 0
    for {set i 0} {$i <100} {incr i} {
        for {set k 0} {$k <3} {incr k} {
            if {[has_feature "ROTATION"]} {
                if {[has_feature "VERLET_STEP4_VELOCITY"]} {
                    set tolerance_o 1.25E-4
                } else {
                    set tolerance_o 1.0E-2
                }
                set do0 [expr abs([lindex [part 0 print omega_body] $k] -exp(-$gamma(0)*[setmd time] /[lindex $J $k]))]
                set do1 [expr abs([lindex [part 1 print omega_body] $k] -exp(-$gamma(1)*[setmd time] /[lindex $J $k]))]
                if { $do0 > $tolerance_o || $do1 > $tolerance_o } {
                    error_exit "Friction Deviation in omega too large. $i $k $do0 $do1"
                }
            }
            if {[has_feature "VERLET_STEP4_VELOCITY"]} {
                set tolerance_v 4.5E-5
            } else {
                set tolerance_v 1.0E-2
            }
            set dv0 [expr abs([lindex [part 0 print v] $k] -exp(-$gamma(0)*[setmd time] / $mass))]
            set dv1 [expr abs([lindex [part 1 print v] $k] -exp(-$gamma(1)*[setmd time] / $mass))]
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
    set temp(0) 2.5
    set temp(1) 2.0
    
    for {set j 0} {$j<3} {incr j} {
        for {set k 0} {$k<2} {incr k} {
            set gamma_tran($k,$j) [expr (0.4 + [t_random]) * 10]
            set gamma_rot($k,$j) [expr (0.2 + [t_random]) * 20]
            if {![has_feature "PARTICLE_ANISOTROPY"]} {
                set gamma_tran($k,1) $gamma_tran($k,0)
                set gamma_tran($k,2) $gamma_tran($k,0)
            }
        
            if {![has_feature "ROTATIONAL_INERTIA"]} {
                set gamma_rot($k,1) $gamma_rot($k,0)
                set gamma_rot($k,2) $gamma_rot($k,0)
            }
        }
    }

    set box 10
    setmd box_l $box $box $box
    set kT 1.5
    
    for {set j 0} {$j<3} {incr j} {
        for {set k 0} {$k<2} {incr k} {
            if {$test_case == 2 || $test_case == 3} {
                set halfkT($k) [expr $temp($k)/2.]
            } else {
                set halfkT($k) [expr $kT / 2.0]
            }
            
            if {$test_case == 1 || $test_case == 3} {
                set gamma_tr($k,$j) $gamma_tran($k,$j)
            } else {
                set gamma_tr($k,$j) $gamma(0)
            }
            # translational diffusion
            set D_tr($k,$j) [expr 2 * $halfkT($k) / $gamma_tr($k,$j)]
        }
    }
    
    if {[has_feature "PARTICLE_ANISOTROPY"]} {
        thermostat langevin $kT $gamma(0) $gamma(0) $gamma(0)
    } else {
        thermostat langevin $kT $gamma(0)
    }

    # no need to rebuild Verlet lists, avoid it
    setmd skin 1.0
    setmd time_step 0.008
    set n 200

    # 0.2 is required to avoid infinitely small requirement [time_step < mass / gamma] and [j / gamma_rot]
    set mass [expr (0.2 + [t_random]) * 7]
    set j1 [expr (0.2 + [t_random]) * 7]
    set j2 [expr (0.2 + [t_random]) * 7]
    set j3 [expr (0.2 + [t_random]) * 7]
    
    if {![has_feature "ROTATIONAL_INERTIA"]} {
        set j2 $j1
        set j3 $j1
    }
    
    for {set i 0} {$i<$n} {incr i} {
        for {set k 0} {$k<2} {incr k} {
            if {[has_feature "ROTATIONAL_INERTIA"]} {
                part [expr $i + $k*$n] pos [expr [t_random] *$box] [expr [t_random] * $box] [expr [t_random] * $box] rinertia $j1 $j2 $j3 mass $mass v 0 0 0
            } else {
                part [expr $i + $k*$n] pos [expr [t_random] *$box] [expr [t_random] * $box] [expr [t_random] * $box] rinertia $j1 mass $mass v 0 0 0
            }
            if {[has_feature "ROTATION"]} {
                part [expr $i + $k*$n] omega_body 0 0 0
            }
            if {[has_feature "PARTICLE_ANISOTROPY"]} {
                    if {[has_feature "ROTATIONAL_INERTIA"]} {
                        switch $test_case {
                            1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2) gamma_rot $gamma_rot($k,0) $gamma_rot($k,1) $gamma_rot($k,2)}
                            2 {part [expr $i + $k*$n] temp $temp($k)}
                            3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2) gamma_rot $gamma_rot($k,0) $gamma_rot($k,1) $gamma_rot($k,2) temp $temp($k)}
                        }
                    } else {
                        if {[has_feature "ROTATION"]} {
                            switch $test_case {
                                1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2) gamma_rot $gamma_rot($k,0)}
                                2 {part [expr $i + $k*$n] temp $temp($k)}
                                3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2) gamma_rot $gamma_rot($k,0) temp $temp($k)}
                            }
                        } else {
                            switch $test_case {
                                1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2)}
                                2 {part [expr $i + $k*$n] temp $temp($k)}
                                3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) $gamma_tran($k,1) $gamma_tran($k,2) temp $temp($k)}
                            }
                        }
                    }
            } else {
                    if {[has_feature "ROTATIONAL_INERTIA"]} {
                        switch $test_case {
                            1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) gamma_rot $gamma_rot($k,0) $gamma_rot($k,1) $gamma_rot($k,2)}
                            2 {part [expr $i + $k*$n] temp $temp($k)}
                            3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) gamma_rot $gamma_rot($k,0) $gamma_rot($k,1) $gamma_rot($k,2) temp $temp($k)}
                        }
                    } else {
                        if {[has_feature "ROTATION"]} {
                            switch $test_case {
                                1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) gamma_rot $gamma_rot($k,0)}
                                2 {part [expr $i + $k*$n] temp $temp($k)}
                                3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) gamma_rot $gamma_rot($k,0) temp $temp($k)}
                            }
                        } else {
                            switch $test_case {
                                1 {part [expr $i + $k*$n] gamma $gamma_tran($k,0)}
                                2 {part [expr $i + $k*$n] temp $temp($k)}
                                3 {part [expr $i + $k*$n] gamma $gamma_tran($k,0) temp $temp($k)}
                            }
                        }
                    }
            }
        }
    }

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

    set loops 100
    puts "Thermalizing..."
    set therm_steps 1200
    integrate $therm_steps

    set int_steps 100
    for {set i 0} {$i <$loops} {incr i} {
        integrate $int_steps
        # Get kinetic energy in each degree of freedom for all particles
        for {set p 0} {$p <$n} {incr p} {
            for {set k 0} {$k<2} {incr k} {
                set ind [expr $p + $k*$n]
                set pos [part [expr $p + $k*$n] print pos]
                set v [part [expr $p + $k*$n] print v]
                if {[has_feature "ROTATION"]} {
                    set o [part [expr $p + $k*$n] print omega_body]
                    set ox2($k) [expr $ox2($k) +pow([lindex $o 0],2)]
                    set oy2($k) [expr $oy2($k) +pow([lindex $o 1],2)]
                    set oz2($k) [expr $oz2($k) +pow([lindex $o 2],2)]
                }
                set vx2($k) [expr $vx2($k) +pow([lindex $v 0],2)]
                set vy2($k) [expr $vy2($k) +pow([lindex $v 1],2)]
                set vz2($k) [expr $vz2($k) +pow([lindex $v 2],2)]
                set dx2($k) [expr pow([expr [lindex $pos 0] - [lindex $pos0($ind) 0]], 2)]
                set dy2($k) [expr pow([expr [lindex $pos 1] - [lindex $pos0($ind) 1]], 2)]
                set dz2($k) [expr pow([expr [lindex $pos 2] - [lindex $pos0($ind) 2]], 2)]
                set dt0(0) [expr $mass / $gamma_tr($k,0)]
                set dt0(1) [expr $mass / $gamma_tr($k,1)]
                set dt0(2) [expr $mass / $gamma_tr($k,2)]
                set dt [expr ($int_steps * ($i + 1) + $therm_steps) * [setmd time_step]]
                # translational diffusion variance: after a closed-form integration of the Langevin EOM
                set sigma2x_tr($k) [expr $D_tr($k,0) * (2 * $dt + $dt0(0) * (-3 + 4 * exp(-$dt / $dt0(0)) - exp(-2 * $dt / $dt0(0))))]
                set sigma2y_tr($k) [expr $D_tr($k,1) * (2 * $dt + $dt0(1) * (-3 + 4 * exp(-$dt / $dt0(1)) - exp(-2 * $dt / $dt0(1))))]
                set sigma2z_tr($k) [expr $D_tr($k,2) * (2 * $dt + $dt0(2) * (-3 + 4 * exp(-$dt / $dt0(2)) - exp(-2 * $dt / $dt0(2))))]
                set sigma2_tr($k) [expr $sigma2x_tr($k) + $sigma2y_tr($k) + $sigma2z_tr($k)]
                set dr($k) [expr $dr($k) + (($dx2($k)+$dy2($k)+$dz2($k)) - $sigma2_tr($k)) / $sigma2_tr($k)]
            }
        }
    }

    set tolerance 0.15

    for {set k 0} {$k<2} {incr k} {
        set Evx($k) [expr 0.5 * $mass *$vx2($k)/$n/$loops]
        set Evy($k) [expr 0.5 * $mass *$vy2($k)/$n/$loops]
        set Evz($k) [expr 0.5 * $mass *$vz2($k)/$n/$loops]

        if {[has_feature "ROTATION"]} {
            set Eox($k) [expr 0.5 * $j1 *$ox2($k)/$n/$loops]
            set Eoy($k) [expr 0.5 * $j2 *$oy2($k)/$n/$loops]
            set Eoz($k) [expr 0.5 * $j3 *$oz2($k)/$n/$loops]
            set do($k) [expr 1./3. *($Eox($k) +$Eoy($k) +$Eoz($k))/$halfkT($k)-1.]
            set dox($k) [expr ($Eox($k))/$halfkT($k)-1.]
            set doy($k) [expr ($Eoy($k))/$halfkT($k)-1.]
            set doz($k) [expr ($Eoz($k))/$halfkT($k)-1.]
        }
        
        set dv($k) [expr 1./3. *($Evx($k) +$Evy($k) +$Evz($k))/$halfkT($k)-1.]
        set dr($k) [expr $dr($k)/$n/$loops]
        
        
        puts "\n"
        puts "1/2 kT = $halfkT($k)"
        puts "translation: $Evx($k) $Evy($k) $Evz($k)"
        if {[has_feature "ROTATION"]} {
            puts "rotation: $Eox($k) $Eoy($k) $Eoz($k)"
            puts "Deviation in rotational energy: $do($k)"
            puts "Deviation in rotational energy per degrees of freedom: $dox($k) $doy($k) $doz($k)"
            if { abs($do($k)) > $tolerance } {
                puts "Moment of inertia principal components: $j1 $j2 $j3"
                error "Relative deviation in rotational energy too large: $do($k)"
            }
        }

        puts "Deviation in translational energy: $dv($k)"
        puts "Deviation in translational diffusion: $dr($k) $mass $gamma_tr($k,0) $gamma_tr($k,1) $gamma_tr($k,2)"

        if { abs($dv($k)) > $tolerance } {
           error "Relative deviation in translational energy too large: $dv($k)"
        }
        if { abs($dr($k)) > $tolerance } {
           error "Relative deviation in translational diffusion too large: $dr($k) for parameters: mass=$mass gamma_tr=$gamma_tr($k,0) $gamma_tr($k,1) $gamma_tr($k,2)"
        }
       
       if {[has_feature "ROTATION"]} {
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
}

# the actual testing
for {set i 0} {$i < 4} {incr i} {
    test_mass-and-rinertia_per_particle $i
}

exit 0
