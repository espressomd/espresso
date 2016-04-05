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
source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "EWALD_GPU"
# for more than 2 nodes, the cutoff will become too large
require_max_nodes_per_side 2

puts "-------------------------------------------"
puts "- Testcase ewaldgpu.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------"


set epsilon_energy 1e-2
set epsilon_force 1e-4
thermostat off
setmd time_step 0.01
setmd skin 0.05

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

if { [catch {
    read_data "ewaldgpu_system.data"
    set target_energy -5168.685115780581

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    for {set testcase 0} {1} {incr testcase} {
        set lb 1.0
        if {$testcase == 0} {
            set testname "Original ewaldgpu testcase"
            set should_succeed 1
            inter coulomb $lb ewaldgpu 11.288055777549744 8 0.2998850205540657
            set intercoulomb {coulomb 1.0 ewaldgpu 11.288055777549744 8 8 8 0.2998850205540657}
        } elseif {$testcase == 1} {
            set testname "Parameter changes"
            set should_succeed 0
            inter coulomb $lb ewaldgpu 11.0 8 0.3
            inter coulomb $lb ewaldgpu 12.0 9 0.5
            set intercoulomb {coulomb 1.0 ewaldgpu 12.0 9 9 9 0.5}
        } elseif {$testcase == 2} {
            set testname "Tuning parameters"
            set should_succeed 1
            inter coulomb $lb ewaldgpu tune accuracy 1e-4 K_max 15 precision 0.00001
            set intercoulomb {coulomb 1.0 ewaldgpu tune accuracy 1e-4 K_max 15 precision 0.00001}
        } elseif {$testcase == 3} {
            set testname "Disabling ewaldgpu"
            set should_succeed 1
            set lb 0.0
            inter coulomb $lb
            set intercoulomb "coulomb 0.0"
        } elseif {$testcase == 4} {
            set testname "Turning ewaldgpu back on"
            set should_succeed 1
            inter coulomb $lb ewaldgpu 11.288055777549744 8 0.2998850205540657
            set intercoulomb {coulomb 1.0 ewaldgpu 11.288055777549744 8 8 8 0.2998850205540657}
        } elseif {$testcase == 5} {
            set testname "Changing Bjerrum length"
            set should_succeed 1
            set lb 2.0
            inter coulomb $lb ewaldgpu 11.288055777549744 8 0.2998850205540657
            set intercoulomb {coulomb 2.0 ewaldgpu 11.288055777549744 8 8 8 0.	}
        } else {
            break
        }

        puts "+++ Running $testname +++"

        # to ensure force recalculation
        integrate 0 recalc_forces
        
        # Only applicable for mmm1dgpu because different tuning results in ewaldgpu possible
					#set tuning_result [inter coulomb]
					#if { [string index $tuning_result 0] == "\{" } { 
					#    set tuning_result [string range $tuning_result 1 end-1]
					#}
					#set tuning_result [split $tuning_result]
					#if { $tuning_result != $intercoulomb } {
					#    if {[lindex $tuning_result 0] != [lindex $intercoulomb 0] || [lindex $tuning_result 2] != [lindex $intercoulomb 2]} {
					#        error "Expected tuning result $intercoulomb, got $tuning_result"
					#    }
					#    foreach a [list 1	] {
					#        if { [lindex $intercoulomb $a] == "ignore" } { continue }
					#        if {abs([lindex $tuning_result $a] - [lindex $intercoulomb $a]) > $epsilon} {
					#            error "Expected tuning result $intercoulomb, got $tuning_result"
					#        }
					#    }
					#}

        set energy [lindex [lindex [analyze energy] 0] 1]
        set dE [expr abs($energy - [expr $target_energy*$lb])]

        set maxdx 0
        set maxpx 0
        set maxdy 0
        set maxpy 0
        set maxdz 0
        set maxpz 0
        for { set i 0 } { $i <= [setmd max_part] } { incr i } {
        	set resF [part $i pr f]
        	#set tgtF $F($i)
          set tgtF [vecscale $lb $F($i)]
        	set dx [expr abs([lindex $resF 0] - [lindex $tgtF 0])]
        	set dy [expr abs([lindex $resF 1] - [lindex $tgtF 1])]
        	set dz [expr abs([lindex $resF 2] - [lindex $tgtF 2])]

        	if { $dx > $maxdx} {
        	    set maxdx $dx
        	    set maxpx $i
        	}
        	if { $dy > $maxdy} {
        	    set maxdy $dy
        	    set maxpy $i
        	}
        	if { $dz > $maxdz} {
        	    set maxdz $dz
        	    set maxpz $i
        	}
        }

        set failed 0
        if { $should_succeed == 1} {
            puts "maximal force deviation in x $maxdx for particle $maxpx, in y $maxdy for particle $maxpy, in z $maxdz for particle $maxpz"
            puts "energy deviation $dE"
        }
        if { $maxdx > $epsilon_force || $maxdy > $epsilon_force || $maxdz > $epsilon_force || $dE > $epsilon_energy } {
            if { $should_succeed == 1 } {
                if { $maxdx > $epsilon_force} {puts "force of particle $maxpx: [part $maxpx pr f] != [vecscale $lb $F($maxpx)]"}
                if { $maxdy > $epsilon_force} {puts "force of particle $maxpy: [part $maxpy pr f] != [vecscale $lb $F($maxpy)]"}
                if { $maxdz > $epsilon_force} {puts "force of particle $maxpz: [part $maxpz pr f] != [vecscale $lb $F($maxpz)]"}
                if { $dE > $epsilon_energy} {puts "Energy: $energy != [expr $lb*$target_energy]"}
                error "error too large on $testname"
            } elseif {($maxdx > $epsilon_force || $maxdy > $epsilon_force || $maxdz > $epsilon_force) && $dE > $epsilon_energy} {
                # if we expect it to fail, both forces and energies have to fail
                set failed 1
            }
        }
        if { $should_succeed == 0 && $failed == 0 } {
            error "$testname should have failed"
        }
        puts [inter coulomb] 
        puts "OK"
    }
} res ] } {
    error_exit $res
}

exit 0
