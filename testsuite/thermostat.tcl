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

set epsilon 3e-2
thermostat langevin 1.0 1.0
setmd time_step 0.01
setmd skin 0.5
set n_part 100
set maxstep 100
set intstep [expr round(50*$maxstep/100)]

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

    set eng0    [analyze energy kin]
    set temp0   [expr $eng0/$n_part/($deg_free/2.)]
    set curtemp1 0

    for {set i 0} { $i < 100} { incr i } {
    integrate $intstep
	puts -nonewline " $i percent done\r"
	flush stdout
	set toteng [analyze energy total]
	set cureng [analyze energy kin] 
	set curtemp [expr $cureng/$n_part/($deg_free/2.)] 

	if { [expr abs($toteng - $cureng)] > $epsilon } {
	    error "system has unwanted energy contributions"
	}
	set curtemp1 [expr $curtemp1 + $curtemp]
    }
    set curtemp1 [expr $curtemp1/$maxstep]
    # here you can create a new snapshot
    # write_data $filename


    set rel_temp_error [expr abs(( [setmd temp] - $curtemp1)/$curtemp1)]
    puts "thermostat temperature:          [setmd temp]"
    puts "measured temperature:            $curtemp1"
    puts "relative temperature deviations: $rel_temp_error"
    if { $rel_temp_error > $epsilon } {
	error "relative temperature error too large"
    }
  
} res ] } {
    error_exit $res
}

exit 0
