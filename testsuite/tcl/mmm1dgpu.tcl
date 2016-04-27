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

require_feature "MMM1D_GPU"
require_feature "PARTIAL_PERIODIC"

puts "-------------------------------------------"
puts "- Testcase mmm1dgpu.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------"

cellsystem nsquare

set epsilon 1e-4
thermostat off
setmd time_step 0.01
setmd skin 0.05

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

if { [catch {
    read_data "mmm1d_system.data"
    set target_energy -7.156365298205383
    setmd periodic 0 0 1

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    for {set testcase 0} {1} {incr testcase} {
        set lb 1.0
        if {$testcase == 0} {
            set testname "Original MMM1D testcase"
            set should_succeed 1
            inter coulomb $lb mmm1dgpu 6.0 3 0.0001
            set intercoulomb {coulomb 1.0 mmm1dgpu 6.0 3 1e-4}
        } elseif {$testcase == 1} {
            set testname "Parameter changes"
            set should_succeed 0
            inter coulomb $lb mmm1dgpu 1.0 1 0.1
            inter coulomb $lb mmm1dgpu 2.0 2 0.2
            set intercoulomb {coulomb 1.0 mmm1dgpu 2.0 2 0.2}
        } elseif {$testcase == 2} {
            set testname "Tuning both parameters"
            set should_succeed 1
            inter coulomb $lb mmm1dgpu tune 1e-4
            set intercoulomb {coulomb 1.0 mmm1dgpu ignore ignore 1e-4}
        } elseif {$testcase == 3} {
            set testname "Tuning only bessel cutoff"
            set should_succeed 1
            inter coulomb $lb mmm1dgpu 6.0 1e-4
            set intercoulomb {coulomb 1.0 mmm1dgpu 6.0 3 1e-4}
        } elseif {$testcase == 4} {
            set testname "Disabling MMM1DGPU"
            set should_succeed 1
            set lb 0.0
            inter coulomb $lb
            set intercoulomb "coulomb 0.0"
        } elseif {$testcase == 5} {
            set testname "Turning MMM1D back on"
            set should_succeed 1
            inter coulomb $lb mmm1dgpu 6.0 3 0.0001
            set intercoulomb {coulomb 1.0 mmm1dgpu 6.0 3 1e-4}
        } elseif {$testcase == 6} {
            set testname "Changing Bjerrum length"
            set should_succeed 1
            set lb 2.0
            inter coulomb $lb mmm1dgpu 6.0 3 0.0001
            set intercoulomb {coulomb 2.0 mmm1dgpu 6.0 3 1e-4}
        } else {
            break
        }

        puts "+++ Running $testname +++"

        # to ensure force recalculation
        integrate 0 recalc_forces
        
        set tuning_result [inter coulomb]
        if { [string index $tuning_result 0] == "\{" } { 
            set tuning_result [string range $tuning_result 1 end-1]
        }
        set tuning_result [split $tuning_result]
        if { $tuning_result != $intercoulomb } {
            if {[lindex $tuning_result 0] != [lindex $intercoulomb 0] || [lindex $tuning_result 2] != [lindex $intercoulomb 2]} {
                error "Expected tuning result $intercoulomb, got $tuning_result"
            }
            foreach a [list 1 3 4 5] {
                if { [lindex $intercoulomb $a] == "ignore" } { continue }
                if {abs([lindex $tuning_result $a] - [lindex $intercoulomb $a]) > $epsilon} {
                    error "Expected tuning result $intercoulomb, got $tuning_result"
                }
            }
        }

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
        }
        if { $maxdx > $epsilon || $maxdy > $epsilon || $maxdz > $epsilon || $dE > $epsilon } {
            if { $should_succeed == 1 } {
                if { $maxdx > $epsilon} {puts "force of particle $maxpx: [part $maxpx pr f] != [vecscale $lb $F($maxpx)]"}
                if { $maxdy > $epsilon} {puts "force of particle $maxpy: [part $maxpy pr f] != [vecscale $lb $F($maxpy)]"}
                if { $maxdz > $epsilon} {puts "force of particle $maxpz: [part $maxpz pr f] != [vecscale $lb $F($maxpz)]"}
                if { $dE > $epsilon} {puts "Energy: $energy != [expr $lb*$target_energy]"}
                error "error too large on $testname"
            } elseif {($maxdx > $epsilon || $maxdy > $epsilon || $maxdz > $epsilon) && $dE > $epsilon} {
                # if we expect it to fail, both forces and energies have to fail
                set failed 1
            }
        }
        if { $should_succeed == 0 && $failed == 0 } {
            error "$testname should have failed"
        }
        puts "OK"
    }

    set testname "Analytical result for R=0, z=1"
    puts "+++ Running $testname +++"
    setmd box_l 10 10 10
    part deleteall
    part 0 pos 0 0 0 q 1
    part 1 pos 0 0 1 q 1
    inter coulomb $lb mmm1dgpu 6.0 3 0.0001
    integrate 0 recalc_forces
    set energy [lindex [lindex [analyze energy] 0] 1]
    if { [lindex [part 0 print f] 0] != 0.0 || [lindex [part 0 print f] 1] != 0.0 
        || [expr abs([lindex [part 0 print f] 2] + 0.99510759)] > $epsilon
        || abs([lindex [part 0 print f] 2] +[lindex [part 1 print f] 2]) > $epsilon
        || [expr abs($energy - 1.00242505606)] > $epsilon
        } {
            error "$testname failed"
    }
    puts "OK"

    set testname "Comparing near/far formula"
    puts "+++ Running $testname +++"
    for {set i 0} {$i < 2} {incr i} {
        part deleteall
        setmd box_l 10 10 10
        if {$i == 0} {
            part 0 pos 0 0 0 q 1
            part 1 pos 4.0 5.0 0.8 q 1
            set switch_near 9.0
            set switch_far 3.0
            set target_near {-0.021209673234069693 -0.026512091542587116 -0.0013773300981728565}
            set target_far {-0.021209673234069689 -0.026512091542587112 -0.0013773300981728574}
        } elseif {$i == 1} {
            part 0 pos 2.4545 1.89461 -4.1594 q -1
            part 1 pos 2.0592 5.3135 4.6326 q -1
            set switch_near 4.0
            set switch_far 3.0
            set target_near {0.0090178535955331984 -0.077994053830590729 0.020480946151392336}
            set target_far {0.0090178535955331984 -0.077994053830590729 0.020480946151392329}
        }

        set rxy2 [expr pow([lindex [part 0 print folded_position] 0]-[lindex [part 1 print folded_position] 0],2) + pow([lindex [part 0 print folded_position] 1]-[lindex [part 1 print folded_position] 1],2)]
        if {$rxy2 > [expr pow($switch_near,2)]} {
            error "Not using near formula!"
        } elseif {$rxy2 < [expr pow($switch_far,2)]} {
            error "Not using far formula!"
        }

        # use far formula
        inter coulomb $lb mmm1dgpu $switch_far 1e-4
        integrate 0 recalc_forces
        set force_far [part 0 print f]

        # use near formula
        inter coulomb $lb mmm1dgpu $switch_near 1e-4
        integrate 0 recalc_forces
        set force_near [part 0 print f]
        
        for { set j 0} {$j < 3} {incr j} {
	        if { [expr abs([lindex $target_far $j] - [lindex $force_far $j])] > $epsilon } {
	        	puts "$target_far $force_far"
	       	}
	        if { [expr abs([lindex $target_near $j] - [lindex $force_near $j])] > $epsilon } {
	        	error "Near result incorrect"
	       	}
	        if { [expr abs([lindex $force_far $j] - [lindex $force_near $j])] > $epsilon } {
	        	error "Near and far formulas do not match"
	       	}
        }
        puts "OK $i"
    }


} res ] } {
    error_exit $res
}

exit 0
