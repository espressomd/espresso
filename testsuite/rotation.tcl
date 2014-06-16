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
# 
source "tests_common.tcl"

require_feature "ROTATION"
require_feature "GAY_BERNE"

puts "----------------------------------------------"
puts "- Testcase rotation.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

set epsilon 5e-4
thermostat off

setmd time_step 0.001
proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

if { [catch {
    read_data "gb_system.data"

    setmd skin 0.5

    inter 0 0 gay-berne 1.0 1.0 4.0 3.0 5.0 2.0 1.0
    
    set GBeng_0 [expr [analyze energy gb 0 0]]
    set toteng_0 [analyze energy total]
    if { [expr abs($toteng_0 - $GBeng_0)] > $epsilon } {
	error "system has unwanted energy contribution, i.e. U_GB != U_total"
    }
    puts "energy before integration: [analyze energy]"
  
    integrate 50

    puts "energy after integration: [analyze energy]"

    # check the conservation of the total energy
    set toteng [analyze energy total]
    set rel_eng_error [expr abs(($toteng_0 - $toteng)/$toteng)]
    puts "total energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }
    
    # check new GB energy against expected value
    set GB_expected -2971.72
    set GBeng [expr [analyze energy gb 0 0]]
    set rel_eng_error [expr abs(($GBeng - $GB_expected)/$toteng)]
    puts "   GB energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }

} res ] } {
    error_exit $res
}

exit 0
