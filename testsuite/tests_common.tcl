#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 

proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

proc require_feature {feature {off ""}} {
    global errf
    if {($off == ""    && ! [regexp $feature [code_info]]) ||
	($off == "off" &&   [regexp $feature [code_info]])} {
	if {$off == ""} {
	    puts stderr "wanted feature not compiled in: $feature"
	} {
	    puts stderr "unwanted feature compiled in: $feature"
	}
	exec rm -f $errf
	exit -42
    }
}

proc require_max_nodes_per_side {n} {
    global errf
    foreach s [setmd node_grid] {
	if {$s > $n} {
	    puts stderr "cannot run on [setmd n_nodes] processors,"
	    puts stderr "since max number of nodes per side is $n,"
	    puts stderr "but node grid is [setmd node_grid]"
	    exec rm -f $errf
	    exit -42
	}
    }
}