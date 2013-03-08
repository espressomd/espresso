# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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

set TEST_FAIL 66
set TEST_IGNORE 42
set TEST_OK 0

proc error_exit {{error 0}} {
    global TEST_FAIL
    if { $error != 0 } { puts $error }
    exit $TEST_FAIL
}

proc ignore_exit {{error 0}} {
    global TEST_IGNORE
    if { $error != 0 } { puts $error }
    exit $TEST_IGNORE
}

proc ok_exit {} {
    global TEST_OK
    exit $TEST_OK
}

proc test_catch {script} {
    if { [catch [uplevel 1 $script] res] } then {
	exit $TEST_FAIL
    }
}

proc has_feature {feature {off ""}} {
    if {($off == ""    && ! [regexp "{ $feature }" [code_info]]) ||
	($off == "off" &&   [regexp "{ $feature }" [code_info]])} {
	return 0
    } else {
	return 1
    }
}

proc require_feature {feature {off ""}} {
    if {($off == ""    && ! [regexp "{ $feature }" [code_info]]) ||
	($off == "off" &&   [regexp "{ $feature }" [code_info]])} {
	if {$off == ""} then {
	    ignore_exit "Feature $feature is not activated."
	} else {
	    ignore_exit "Unwanted feature $feature is activated."
	}
    }
}

proc require_max_nodes_per_side {n} {
    foreach s [setmd node_grid] {
	if {$s > $n} {
	    ignore_exit "Testcase cannot run on [setmd n_nodes] processors, \n\tsince max number of nodes per side is $n,\n\tbut node grid is [setmd node_grid]"
	}
    }
}
