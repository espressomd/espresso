# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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


# check the charge-charge P3M  algorithm
source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "EWALD_GPU"

#cellsystem domain_decomposition
#cellsystem nsquare

puts "---------------------------------------------------------------"
puts "- Testcase ewaldgpu.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

set epsilon 1e-3
thermostat off
setmd skin 0.05

set int_steps 10
setmd time_step 0.00001

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos q f}
    close $f
}

if { [catch {
    puts "START"
    puts "List:[cuda list]"
		puts "Device:[cuda getdevice]"
    read_data "ewaldgpu_system.data"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
		set F($i) [part $i pr f]
    }

# INTEGRATION
#inter coulomb 1.0 ewaldgpu tunealpha 11.431538115739823 6 0.0001
inter coulomb 1.0 ewaldgpu tune accuracy 1e-4 K_max 15 precision 0.00001
#inter coulomb 1.0 p3m tune accuracy 1e-4


puts "TUNED"

set start_time [expr 1.0*[clock clicks -milliseconds]]
integrate $int_steps recalc_forces
set end_time [expr 1.0*[clock clicks -milliseconds]]
puts "TIME: [expr $end_time-$start_time] millisec"

write_data "ewaldgpu_system.data.out"

puts "E=[analyze energy kinetic]"
puts "E=[analyze energy coulomb]"
puts "E=[analyze energy]"
puts [inter coulomb]

#end this part of the p3m-checks by cleaning the system .... 
part deleteall
inter coulomb 0.0

} res ] } {
    error_exit $res
}

exit 0