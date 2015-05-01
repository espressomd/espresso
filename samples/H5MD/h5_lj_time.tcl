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
setmd time_step 0.05
setmd skin 0
set time_steps 10

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

inter 0 0 lennard-jones 1.0 1.0 1.12246

#create h5 file
h5mdfile H5Fcreate "h5mdfile.h5" 

h5mdfile H5Gcreate2 "group1"
h5mdfile H5Screate_simple type double dims 1 600 3
h5mdfile H5Pset_chunk dims 1 600 3
h5mdfile H5Dcreate2 "/group1/dset1"
h5mdfile H5Screate_simple type double dims 1 600 3
h5mdfile H5Pset_chunk dims 1 600 3
h5mdfile H5Dcreate2 "/group1/dset2"

h5mdfile H5Gcreate2 "group2"
h5mdfile H5Screate_simple type double dims 1 1
h5mdfile H5Pset_chunk dims 1 1
h5mdfile H5Dcreate2 "/group2/dset3"

#Start value
integrate 0
h5mdfile H5Dopen2 "/group1/dset1"
for { set j 0 } { $j <= [setmd max_part] } { incr j } {
    set pos [part $j pr p]
    h5mdfile H5_write_value value [lindex $pos 0] index 0 $j 0
    h5mdfile H5_write_value value [lindex $pos 1] index 0 $j 1
    h5mdfile H5_write_value value [lindex $pos 2] index 0 $j 2
}
puts kat
h5mdfile H5Dwrite

h5mdfile H5Dopen2 "/group1/dset2"
for { set j 0 } { $j <= [setmd max_part] } { incr j } {
    set force [part $j pr f]
    h5mdfile H5_write_value value [expr [lindex $force 0]] index 0 $j 0
    h5mdfile H5_write_value value [expr [lindex $force 1]] index 0 $j 1
    h5mdfile H5_write_value value [expr [lindex $force 2]] index 0 $j 2
}
h5mdfile H5Dwrite



#integrate
for { set i 0 } { $i < $time_steps } { incr i } {
    integrate 10
    h5mdfile H5Dopen2 "/group1/dset1"
    h5mdfile H5Dextend dims [expr $i+2] 600 3
    h5mdfile H5Sselect_hyperslab offset [expr $i+1] 0 0
    h5mdfile H5Screate_simple type double dims 1 600 3
    for { set j 0 } { $j <= [setmd max_part] } { incr j } {
	set pos [part $j pr p]
	h5mdfile H5_write_value value [lindex $pos 0] index 0 $j 0
	h5mdfile H5_write_value value [lindex $pos 1] index 0 $j 1
	h5mdfile H5_write_value value [lindex $pos 2] index 0 $j 2
    }
    h5mdfile H5Dwrite
    
    h5mdfile H5Dopen2 "/group1/dset2"
    h5mdfile H5Dextend dims [expr $i+2] 600 3
    h5mdfile H5Sselect_hyperslab offset [expr $i+1] 0 0
    h5mdfile H5Screate_simple type double dims 1 600 3
    for { set j 0 } { $j <= [setmd max_part] } { incr j } {
	set force [part $j pr f]
	h5mdfile H5_write_value value [expr [lindex $force 0]] index 0 $j 0
	h5mdfile H5_write_value value [expr [lindex $force 1]] index 0 $j 1
	h5mdfile H5_write_value value [expr [lindex $force 2]] index 0 $j 2
    }
    h5mdfile H5Dwrite
    
    h5mdfile H5Dopen2 "/group2/dset3"
    h5mdfile H5Dextend dims [expr $i+2] 1
    h5mdfile H5Sselect_hyperslab offset [expr $i+1] 0
    h5mdfile H5Screate_simple type double dims 1 1
    h5mdfile H5_write_value value $i index 0 0
    h5mdfile H5Dwrite
}

write_data "h5_lj_time.out"

#close file
h5mdfile H5Dclose
h5mdfile H5Sclose
h5mdfile H5Fclose


puts "END"

} res ] } {
    error_exit $res
}

exit 0


