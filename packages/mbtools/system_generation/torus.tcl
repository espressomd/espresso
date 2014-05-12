# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#   
# ESPResSo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Routines for placing lipids uniformly on a torus
#
#

namespace eval ::mbtools::system_generation {}

namespace eval ::mbtools::system_generation::torus {

    namespace export create_torus
}

proc ::mbtools::system_generation::torus::find_deltafudge { r1 r2 delta_in nlipids } {
    set pi 3.141592
    set deltafudge 0.01
    set fudgechange 0.001
    set ncount 0
    set changecount 0
    while { $ncount != $nlipids && $changecount < 500 } {
	set delta [expr $delta_in + $deltafudge]    
	set mtheta [expr int(floor(2*$pi*$r2/sqrt($delta)))]
	set dtheta [expr 2*$pi*$r2/(1.0*$mtheta)]
	set dphi [expr $delta/$dtheta]
	set ncount 0
	for { set m 0 } { $m < $mtheta } { incr m } {
	    set theta [expr 2*$pi*$m/($mtheta*1.0)]
	    set mphi [expr int(floor(2*$pi*($r1-$r2*cos($theta))/$dphi))]
	    for { set n 0 } { $n < $mphi } { incr n } {
#		set phi [expr 2*$pi*$n/(1.0*$mphi)]
		incr ncount
	    }
	}

	if { $ncount > $nlipids && $fudgechange < 0 } {
	    set fudgechange [expr -$fudgechange/(2.0) ]
	    incr changecount
	}
	if { $ncount < $nlipids && $fudgechange > 0 } {
	    set fudgechange [expr -$fudgechange/(2.0) ]
	    incr changecount
	}

	set deltafudge [expr $deltafudge + $fudgechange ]

#	puts "$ncount $fudgechange $deltafudge"
    }
    set deltafudge [expr $deltafudge - $fudgechange ]
    return $deltafudge 
 
}

proc ::mbtools::system_generation::torus::create_torus { args } {

    mmsg::send [namespace current] "placing lipids on a torus "
    # ---- Process Command Line Args --------------------------- #
    
    set options {
	{c.arg  { 0.0 0.0 0.0 } "location of the center of the torus relative to the box center" }
	{initarea.arg      1.29    "the starting value for area per mol" }
	{shuffle "shuffle the topology before placing molecules "}
	{bondl.arg     1.0   "bond length between atoms"  }
	{ratio.arg 1.414213 "ratio between inner and outer torus radii"}
    }
    set usage "Usage: create_torus \[c:initarea:shuffle:bondl:ratio\] ]: "
    array set params [::cmdline::getoptions args $options $usage]

    global ::mbtools::system_generation::moltypeskey
    global ::mbtools::system_generation::boxl
    global ::mbtools::system_generation::topology

    if { $params(shuffle) } {
	set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }

    set area_lipid $params(initarea)
    set ratio $params(ratio)

    #Construct the center
    set center [list 0.0 0.0 0.0 ]
    lset center 0 [expr [lindex $params(c) 0] + [lindex $boxl 0]/(2.0)]
    lset center 1 [expr [lindex $params(c) 1] + [lindex $boxl 1]/(2.0)]
    lset center 2 [expr [lindex $params(c) 2] + [lindex $boxl 2]/(2.0)]


    mmsg::send [namespace current] "creating torus with ratio: $ratio at $center"


    # First work out how many mol types there are and construct a list
    # of their lengths
    
    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]


    # Work out the average molecule size
    set sum 0
    foreach mol $topology {
	set sum [expr $sum + [llength $mol] - 1]
    }
    set avmolsize [expr $sum/(1.0*$nmols)]

    set pi 3.141592
    set del [expr $avmolsize/2.0]

    set rad2_tot [expr sqrt($nmols*$area_lipid/(8.0*$pi*$pi*$ratio)) ]
    set rad1 [expr $rad2_tot*$ratio]

    mmsg::send [namespace current] "r1: $rad1 r2: $rad2_tot"

    set nmols_outer [expr 4*$pi*$pi*$rad1*($rad2_tot+$del)/$area_lipid]
    set nmols_inner [expr 4*$pi*$pi*$rad1*($rad2_tot-$del)/$area_lipid]

    set nmols_outer [expr int(($nmols_outer))]
    set nmols_inner [expr $nmols - $nmols_outer]
     mmsg::send [namespace current]"N+: $nmols_outer N-:$nmols_inner"


    # Determine r1 and r2 for the outer layer
    set rad2 [expr $rad2_tot + $del ]
    mmsg::send [namespace current] "r1: $rad1 r2: $rad2"
 
    set deltafudge [ find_deltafudge $rad1 $rad2 $area_lipid $nmols_outer ]
    set delta [expr $area_lipid + $deltafudge]
    
    mmsg::send [namespace current] "Old Area_lipid: $area_lipid New Area_lipid: $delta"

    # Now go through and place all lipids in the outer layer
    set mtheta [expr int((2*$pi*$rad2/sqrt($delta)))]
    set dtheta [expr 2*$pi*$rad2/(1.0*$mtheta)]
    set dphi [expr $delta/$dtheta]
    set molnum 0
    set outercount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr 2*$pi*$m/($mtheta*1.0)]
	set mphi [expr int((2*$pi*($rad1-$rad2*cos($theta))/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    lappend tailpos [expr ($rad1 - $rad2*cos($theta))*cos($phi) + [lindex $center 0]]
	    lappend tailpos [expr ($rad1 - $rad2*cos($theta))*sin($phi) + [lindex $center 1]]
	    lappend tailpos [expr $rad2*sin($theta) + [lindex $center 2]]

	    # Find the normal vector
	    lappend orient [expr -cos($theta)*cos($phi)]
	    lappend orient [expr -cos($theta)*sin($phi)]
	    lappend orient [expr sin($theta)]

	    set orient [::mbtools::utils::normalize $orient]
	    for { set x 0 } { $x < 3 } { incr x } {
		lset tailpos $x [expr [lindex $tailpos $x] - ($del-0.5)*[lindex $orient $x]*$params(bondl) ]
	    }



	    # Place the molecule
	    set mol [lindex $topology $molnum]

	    ::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)

	    unset orient
	    unset tailpos

	    incr molnum
	    incr outercount
	}

    }
    mmsg::send [namespace current] "Placed $outercount lipids in outer layer"

    #------------------ Inner ---------#
    # Determine r1 and r2 for the inner layer
    set rad2 [expr $rad2_tot - $del ]
    mmsg::send [namespace current] "r1: $rad1 r2: $rad2"
 
    set deltafudge [ find_deltafudge $rad1 $rad2 $area_lipid $nmols_inner ]
    set delta [expr $area_lipid + $deltafudge]

    # Now go through and place all lipids in the outer layer
    set mtheta [expr int((2*$pi*$rad2/sqrt($delta)))]
    set dtheta [expr 2*$pi*$rad2/(1.0*$mtheta)]
    set dphi [expr $delta/$dtheta]
    set innercount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr 2*$pi*$m/($mtheta*1.0)]
	set mphi [expr int((2*$pi*($rad1-$rad2*cos($theta))/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    lappend tailpos [expr ($rad1 - $rad2*cos($theta))*cos($phi) + [lindex $center 0]]
	    lappend tailpos [expr ($rad1 - $rad2*cos($theta))*sin($phi) + [lindex $center 1]]
	    lappend tailpos [expr $rad2*sin($theta) + [lindex $center 2]]
	    
	    # Find the normal vector pointing inwards
	    lappend orient [expr cos($theta)*cos($phi)]
	    lappend orient [expr cos($theta)*sin($phi)]
	    lappend orient [expr -sin($theta)]
	    
	    set orient [::mbtools::utils::normalize $orient]
	    for { set x 0 } { $x < 3 } { incr x } {
		lset tailpos $x [expr [lindex $tailpos $x] - ($del-0.5)*[lindex $orient $x]*$params(bondl) ]
	    }

	    # Place the molecule
	    set mol [lindex $topology $molnum]

	    ::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)

	    unset orient
	    unset tailpos

	    incr molnum
	    incr innercount
	    
	}
    }
    mmsg::send [namespace current] "Placed $innercount lipids in inner layer"
    #    exit
    
    mmsg::send [namespace current] "Toroidal vesicle created" 
    flush stdout
    
    return 
}


