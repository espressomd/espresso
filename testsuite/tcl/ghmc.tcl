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

puts "------------------------------------------------"
puts "- Testcase thermostat.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

require_feature THERMOSTAT_IGNORE_NON_VIRTUAL off
require_feature "GHMC"

# we expect to stay in this confidence interval (times stddeviation)
# 20 is a really small chance, but since we measure quite a lot,
# that is still quite a small interval, and the test probably only
# fails if there is really something wrong
set confidence 20
set maxstep 100
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
    blockfile $f write particles {id pos v omega_lab} 
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
    set box_l 10
    setmd box_l $box_l $box_l $box_l
    for {set i 0} {$i < 100} {incr i} {
    	part $i pos [expr rand()*$box_l] [expr rand()*$box_l] [expr rand()*$box_l]
    }

    # make real random draw
    set cmd "t_random seed"
    for {set i 0} {$i < [setmd n_nodes]} { incr i } {
	lappend cmd [expr [pid] + $i] }
    eval $cmd


    thermostat ghmc 1.0 10 0.0
    setmd time_step 0.01
    setmd skin 0.5

    # make sure to give the thermostat a chance to heat up
    integrate 1000

    # checking kinetic energy
    set eng0    [analyze energy kin]
    set temp0   [expr $eng0/[setmd n_part]/($deg_free/2.)]
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
	set curtemp [expr $cureng/[setmd n_part]/($deg_free/2.)] 

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
		set curvel2($c) [expr $curvel2($c) + pow($vx, 2)]
		set curvel4($c) [expr $curvel4($c) + pow($vx, 4)]
	    }
	}
    }

    # here you can create a new snapshot
    # write_data $filename

    puts "NOTE: this is a statistical test, which can fail,"
    puts "even if everything works correctly. However, the"
    puts "chance is really SMALL, so if you see an error"
    puts "here, consider a bug in the thermostat.\n"

    puts "checking velocity distribution"
    puts "output is coordinate: moments 1, 2 and 4"
    puts "expected values are close to 0, 1 and 3"
    set accept1 [expr $confidence*2.*1/$maxstep/[setmd n_part]]
    set accept2 [expr $confidence*2.*30/$maxstep/[setmd n_part]]
    puts "maximally accepted 1st moment deviation [expr sqrt($accept1)]"
    puts "maximally accepted 2nd moment deviation [expr sqrt($accept2)]"

    set err_in_mean 0
    set err_in_var 0
    for {set c 0} {$c < 3} {incr c} {
	set vel1 [expr $curvel1($c)/$maxstep/[setmd n_part]]
	set vel2 [expr $curvel2($c)/$maxstep/[setmd n_part]]
	set vel4 [expr $curvel4($c)/$maxstep/[setmd n_part]]
	puts "${c}:\t$vel1\t$vel2\t$vel4"
	if {pow($vel1, 2) > $accept1} {
	    set err_in_mean [expr pow($vel1, 2)]
	}
	if {pow($vel2 - 1, 2) > $accept2} {
	    set err_in_var [expr pow($vel2 - 1, 2)]
	}
    }
    if {$err_in_mean} {
	error "velocity mean extremely large ($err_in_mean vs. $accept1)"
    }
    if {$err_in_var} {
	error "velocity variance extremely large ($err_in_var vs. $accept2)"
    }
    set samples [expr [setmd n_part]*[degrees_of_freedom]*$maxstep]

    set mean    [expr $curtemp1/$maxstep]
    set var     [expr $curtemp2/$maxstep - $mean*$mean]
    set var_dof [expr $var*[setmd n_part]*$deg_free]

    set epsilon [expr $confidence*2.0/$samples]
    set eps_var [expr 30*$epsilon*[setmd n_part]]

    puts "expected temperature:         [setmd temp]"
    puts "measured temperature:         $mean"
    puts "expected deviations are below [expr sqrt($epsilon)]"
    puts "expected variance per DOF:    2"
    puts "variance per DOF:             $var_dof"
    puts "expected deviations are below [expr sqrt($eps_var)]"

    if {pow($mean - [setmd temp], 2) > $epsilon} {
	error "temperature deviates unusually strongly from [setmd temp]"
    }
  
    # factor 3 from 4th moment
    if {pow($var_dof - 2.0, 2) > $eps_var} {
	error "temperature variance extremely large"
    }
    
} res ] } {
    error_exit $res
}

exit 0
