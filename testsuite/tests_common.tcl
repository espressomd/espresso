# Copyright (C) 2010,2011 The ESPResSo project
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

# 

proc error_exit {error} {
    puts stderr "Error occured: $error"
    exit -666
}

proc has_feature {feature {off ""}} {
    if {($off == ""    && ! [regexp "{ $feature }" [code_info]]) ||
	($off == "off" &&   [regexp "{ $feature }" [code_info]])} {
	return 0;
    } else {
	return 1;
    }
}

proc require_feature {feature {off ""}} {
    if {($off == ""    && ! [regexp "{ $feature }" [code_info]]) ||
	($off == "off" &&   [regexp "{ $feature }" [code_info]])} {
	if {$off == ""} {
	    puts stderr "wanted feature not compiled in: $feature"
	} {
	    puts stderr "unwanted feature compiled in: $feature"
	}
	exit -42
    }
}

proc require_max_nodes_per_side {n} {
    foreach s [setmd node_grid] {
	if {$s > $n} {
	    puts stderr "cannot run on [setmd n_nodes] processors,"
	    puts stderr "since max number of nodes per side is $n,"
	    puts stderr "but node grid is [setmd node_grid]"
	    exit -42
	}
    }
}
