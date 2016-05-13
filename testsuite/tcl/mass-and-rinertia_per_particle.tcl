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
puts "- Testcase mass-and-rinertia_per_particle.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"
puts "------------------------------------------------"

proc test_mass-and-rinertia_per_particle {test_case} {
    set gamma0 1.0
    set gamma1 1.0
    
    # Decelleration
    setmd skin 0
    setmd time_step 0.01
    thermostat langevin 0 $gamma0 
    set J "10 10 1"

    part deleteall
    part 0 pos 0 0 0 rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1 
    part 1 pos 0 0 0 rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] omega_body 1 1 1
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

    for {set i 0} {$i <100} {incr i} {
        for {set k 0} {$k <3} {incr k} {
            if { abs([lindex [part 0 print omega_body] $k] -exp(-$gamma0*$i/10. /[lindex $J $k])) >0.01 ||
                 abs([lindex [part 1 print omega_body] $k] -exp(-$gamma1*$i/10. /[lindex $J $k])) >0.01
            } {
                error_exit "Friction Deviation in omega too large. $i $k"
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
    set gamma(0) 1.0
    set gamma(1) 2.0
    set temp(0) 2.5
    set temp(1) 2.0

    set box 10
    setmd box_l $box $box $box
    set kT 1.5
    set halfkT 0.75
    
    thermostat langevin $kT 1

    # no need to rebuild Verlet lists, avoid it
    setmd skin 1.0
    setmd time_step 0.01

    set n 100
    set mass [expr [t_random] *20]
    set j1 [expr [t_random] * 20]
    set j2 [expr [t_random] * 20]
    set j3 [expr [t_random] * 20]

    for {set i 0} {$i<$n} {incr i} {
        for {set k 0} {$k<2} {incr k} {
            part [expr $i + $k*$n] pos [expr [t_random] *$box] [expr [t_random] * $box] [expr [t_random] * $box] rinertia $j1 $j2 $j3 mass $mass
            switch $test_case {
                1 {part [expr $i + $k*$n] gamma $gamma($k)}
                2 {part [expr $i + $k*$n] temp $temp($k)}
                3 {part [expr $i + $k*$n] gamma $gamma($k) temp $temp($k)}
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
    }


    set loops 100
    puts "Thermalizing..."
    integrate 1000
    puts "Measuring..."

    for {set i 0} {$i <$loops} {incr i} {
        integrate 100
        # Get kinetic energy in each degree of freedom for all particles
        for {set p 0} {$p <$n} {incr p} {
            for {set k 0} {$k<2} {incr k} {
                set v [part [expr $p + $k*$n] print v]
                set o [part [expr $p + $k*$n] print omega_body]
                set ox2($k) [expr $ox2($k) +pow([lindex $o 0],2)]
                set oy2($k) [expr $oy2($k) +pow([lindex $o 1],2)]
                set oz2($k) [expr $oz2($k) +pow([lindex $o 2],2)]
                set vx2($k) [expr $vx2($k) +pow([lindex $v 0],2)]
                set vy2($k) [expr $vy2($k) +pow([lindex $v 1],2)]
                set vz2($k) [expr $vz2($k) +pow([lindex $v 2],2)]
            }
        }
    }

    set tolerance 0.1
    for {set k 0} {$k<2} {incr k} {
        set Evx($k) [expr 0.5 * $mass *$vx2($k)/$n/$loops]
        set Evy($k) [expr 0.5 * $mass *$vy2($k)/$n/$loops]
        set Evz($k) [expr 0.5 * $mass *$vz2($k)/$n/$loops]

        set Eox($k) [expr 0.5 * $j1 *$ox2($k)/$n/$loops]
        set Eoy($k) [expr 0.5 * $j2 *$oy2($k)/$n/$loops]
        set Eoz($k) [expr 0.5 * $j3 *$oz2($k)/$n/$loops]

        if {$test_case == 2 || $test_case ==3} {
            set halfkT [expr $temp($k)/2.]
        }

        set dv($k) [expr 1./3. *($Evx($k) +$Evy($k) +$Evz($k))/$halfkT-1.]
        set do($k) [expr 1./3. *($Eox($k) +$Eoy($k) +$Eoz($k))/$halfkT-1.]
        
        puts "\n"
        puts "1/2 kT = $halfkT"
        puts "translation: $Evx($k) $Evy($k) $Evz($k) rotation: $Eox($k) $Eoy($k) $Eoz($k)"

        puts "Deviation in translational energy: $dv($k)"
        puts "Deviation in rotational energy: $do($k)"

        if { abs($dv($k)) > $tolerance } {
           error "Relative deviation in translational energy too large: $dv($k)"
        }
        if { abs($do($k)) > $tolerance } {
           error "Relative deviation in rotational energy too large: $do($k)"
        }
    }
}

# the actual testing
for {set i 0} {$i < 4} {incr i} {
    test_mass-and-rinertia_per_particle $i
}

exit 0
