# Copyright (C) 2010,2012,2013 The ESPResSo project
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

require_feature "COMFIXED"
require_feature "PARTIAL_PERIODIC"
require_max_nodes_per_side 1

puts "----------------------------------------"
puts "- Testcase comfixed.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

set epsilon 1e-4
thermostat off
setmd time_step 1
setmd skin 0

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    global energy pressure
    set f [open $file "w"]
    set energy [analyze energy total]
    set pressure [analyze pressure total]
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write variable box_l
    blockfile $f write particles {id type pos f}
    close $f
}


if { [catch {
    ############## comfixed-specific part

    setmd box_l 10. 10. 10.
    setmd skin  0.5
    setmd time_step 0.005
    setmd periodic  0 0 0
    thermostat langevin 1.0 1.0

    inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    constraint sphere center 5. 5. 5.  radius 4. type 5
    inter 0 5 lennard-jones 1.0 1.0 1.12246 0.25 0 0
    inter 0 0 comfixed 1

    part 0 pos 5 5 5 type 0
    part 1 pos 6 5 5 type 0
    part 2 pos 6 6 4 type 0
    part 3 pos 4 5 7 type 0

    set com_ini [analyze centermass 0]
    integrate 100
    set com_fin [analyze centermass 0]
    set com_dist [bond_length $com_ini $com_fin]

    if { $com_dist > $epsilon } {
	error "center of mass of particles has moved $com_dist"
    }

} res ] } {
    error_exit $res
}

exit 0
