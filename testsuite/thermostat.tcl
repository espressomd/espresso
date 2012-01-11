# Copyright (C) 2010,2011 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
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

puts "------------------------------------------------"
puts "- Testcase thermostat.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

# we expect to stay in this confidence interval (times stddeviation)
# 3 gives a chance of less than 1 % of failure
set confidence 3
set maxstep 200
set intstep 100

# for checkin unwanted energy contributions
set epsilon 1e-5

set tcl_precision 5

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos v omega} 
    close $f
}
if { [catch {
    if { [regexp "ROTATION" [code_info]] } {
	puts "rotation found, 6 degrees of freedom"
	set deg_free 6
	set filename "thermostat_rot.data"
    } else {
	puts "no rotation found, 3 degrees of freedom"
	set deg_free 3
	set filename "thermostat.data"
    }

    read_data $filename

    # generate some particles
    #for {set i 0} {$i < 1} {incr i} {
    # part $i pos 0 0 0
    #}

    set n_part [setmd n_part]

    # make real random draw
    set cmd "t_random seed"
    for {set i 0} {$i < [setmd n_nodes]} { incr i } {
	lappend cmd [expr [pid] + $i] }
    eval $cmd

    thermostat langevin 1.0 1.0
    setmd time_step 0.01
    setmd skin 0.5

    # make sure to give the thermostat a chance to heat up
    integrate 1000

    # checking kinetic energy
    set eng0    [analyze energy kin]
    set temp0   [expr $eng0/$n_part/($deg_free/2.)]
    set curtemp1 0
    set curtemp2 0

    # checking Maxwell distribution
    # 1, 2, and 4th moment of a single particles velocity
    for {set c 0} {$c < 3} {incr c} {
	set curvel1($c)  0
	set curvel2($c)  0
	set curvel4($c)  0
    }

    for {set i 0} { $i < $maxstep} { incr i } {
	integrate $intstep
	set toteng [analyze energy total]
	set cureng [analyze energy kin] 
	set curtemp [expr $cureng/$n_part/($deg_free/2.)] 

	if { [expr abs($toteng - $cureng)] > $epsilon } {
	    error "system has unwanted energy contributions"
	}
	set curtemp1 [expr $curtemp1 + $curtemp]
	set curtemp2 [expr $curtemp2 + $curtemp*$curtemp]

	for {set p 0} {$p < [setmd n_part]} {incr p} {
	    set v [part $p pr v]
	    for {set c 0} {$c < 3} {incr c} {
		set vx [lindex $v $c]
		set curvel1($c) [expr $curvel1($c) + $vx]
		set curvel2($c) [expr $curvel2($c) + $vx**2]
		set curvel4($c) [expr $curvel4($c) + $vx**4]
	    }
	}
    }
    
    for {set c 0} {$c < 3} {incr c} {
	set vel1 [expr $curvel1($c)/$maxstep/$n_part]
	set vel2 [expr $curvel2($c)/$maxstep/$n_part]
	set vel4 [expr $curvel4($c)/$maxstep/$n_part]
	puts "for the ${c}th coordinate, the 1,2 and 4th velocity moments are $vel1, $vel2 and $vel4"
    }

    set mean [expr $curtemp1/$maxstep]
    set stddev [expr sqrt($curtemp2/$maxstep - $mean*$mean)]

    # here you can create a new snapshot
    # write_data $filename

    set rel_temp_error [expr abs(([setmd temp] - $mean)/[setmd temp])]
    puts "thermostat temperature:          [setmd temp]"
    puts "measured temperature:            $mean"
    puts "fluctuations per DOF:            [expr $stddev*sqrt($n_part*$deg_free)]"
    puts "relative temperature deviation:  $rel_temp_error"

    set epsilon [expr $confidence*sqrt(2)/sqrt($n_part*$deg_free*$maxstep)]
    puts "expected interval of deviation:  $epsilon"

    if { $rel_temp_error > $epsilon } {
	error "relative temperature error too large"
    }
  
} res ] } {
    error_exit $res
}

exit 0
