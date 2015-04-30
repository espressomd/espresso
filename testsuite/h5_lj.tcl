# Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

require_feature "LENNARD_JONES"
if {[has_feature "LEES_EDWARDS"]} {
    require_max_nodes_per_side 2
}


puts "----------------------------------------"
puts "- Testcase lj.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

set epsilon 1e-4
thermostat off
setmd time_step 0.000001
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

puts "START"
read_data "lj_system.data"
for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    set F($i) [part $i pr f]
}

#integrate
inter 0 0 lennard-jones 1.0 1.0 1.12246
integrate 10

#create h5 file
h5mdfile H5Fcreate "h5mdfile.h5" 

#write first dataset
h5mdfile H5Screate_simple double dims 1 600 3
h5mdfile H5Dcreate2 "/dset1"
for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    set pos [part $i pr p]
    h5mdfile H5_write_value value [lindex $pos 0] index 0 $i 0
    h5mdfile H5_write_value value [lindex $pos 1] index 0 $i 1
    h5mdfile H5_write_value value [lindex $pos 2] index 0 $i 2
}
h5mdfile H5Dwrite


#write second dataset
h5mdfile H5Screate_simple double dims 1 600 3
h5mdfile H5Dcreate2 "/dset2"
for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    set force [part $i pr p]
    h5mdfile H5_write_value value [expr 2*[lindex $force 0]] index 0 $i 0
    h5mdfile H5_write_value value [expr 2*[lindex $force 1]] index 0 $i 1
    h5mdfile H5_write_value value [expr 2*[lindex $force 2]] index 0 $i 2
}
h5mdfile H5Dwrite


#close file
h5mdfile H5Dclose
h5mdfile H5Sclose
h5mdfile H5Fclose


puts "END"

} res ] } {
    error_exit $res
}

exit 0
