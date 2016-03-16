# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
# This set of routines constructs a cylindrical vesicle where we attempt
# to place the lipids at a uniform density taking account of the
# cylinder's curvature.
#
# Note: in present version a cylinder with PBCs around one direction is constructed
# Author: Vagelis

namespace eval ::mbtools::system_generation {}

namespace eval ::mbtools::system_generation::cylinder {
    namespace export create_cylinder
}

# ::mbtools::system_generation::find_n_cyl --
#
# Calculate the number of molecules that would be required to coat a
# cylinder of radius r given an area <area> per molecule.
#
proc ::mbtools::system_generation::cylinder::find_n_cyl { r box_l area } {
    set pi 3.141592
    set d [expr sqrt($area)]
    set ellz [lindex $box_l 2]
    set mzeta [expr int(($ellz/$d))]
    set dzeta [expr $ellz/$mzeta]
    set dtheta [expr $area/$dzeta]
    set mtheta [expr int((2.0*$pi*$r/$dtheta))]
    set ncount 0
    for { set m 0 } { $m < $mzeta } { incr m } {
	set ellz_curr [expr ($m+0.5)*$ellz/($mtheta*1.0)]
	for { set n 0 } { $n < $mtheta } { incr n } {
	    set theta [expr 2.0*$pi*$n/$mtheta]
	    incr ncount
	}
    }
    return $ncount
}

# ::mbtools::system_generation::find_fudgearea_cyl --
#
# Attempt to find a value for the area per lipid that allows required
# number of lipids to sit uniformly on a cylinder of the given radius and z-size.
#
proc ::mbtools::system_generation::cylinder::find_fudgearea_cyl { r1  box_l init_area nlipids } {
    set init_step 0.01

    set init_count [find_n_cyl $r1 $box_l $init_area ]
    set init_diff [expr $init_count - $nlipids]
    mmsg::debug [namespace current] "init: A: $init_area R: $r1 diff: $init_diff cnt: [find_n_cyl $r1 $box_l $init_area]"
    # Find initial bracket
    if { $init_diff > 0 } {
	set step [expr $init_step]
	set lower_b $init_area
	set area $init_area
	set diff $init_diff
	while { $diff > 0 } {
	    set area [expr $area + $step]
	    set diff [expr  [find_n_cyl $r1 $box_l $area ] - $nlipids]
	}
	set upper_b $area
    } elseif { $init_diff < 0 } {
	set step [expr -$init_step]
	set upper_b $init_area
	set area $init_area
	set diff $init_diff
	while { $diff < 0 } {
	    set area [expr $area + $step]
	    set diff [expr  [find_n_cyl $r1 $box_l $area ] - $nlipids]
	}
	set lower_b $area
    }

    mmsg::debug [namespace current] "bracket: low: $lower_b [find_n_cyl $r1 $box_l $lower_b] up: $upper_b [find_n_cyl $r1 $box_l $upper_b] "

    # bisect until $diff = 0
    set bisectcount 0
    set maxbisects 50
    while { $diff != 0 && $bisectcount < $maxbisects } {
	set area [expr ($upper_b + $lower_b)/2.0]
	set diff [expr  [find_n_cyl $r1 $box_l $area ]  - $nlipids]
	mmsg::debug [namespace current] "bracket: low: $lower_b [find_n_cyl $r1 $box_l $lower_b] up: $upper_b [find_n_cyl $r1 $box_l $upper_b] new: $area [find_n_cyl $r1 $box_l $area] "
	if { $diff > 0 } {
	    set lower_b $area
	} elseif { $diff < 0 } {
	    set upper_b $area
	}
	incr bisectcount
    }
    if { $bisectcount == $maxbisects } {
	return $lower_b
    } else {
	return $area
    }
}

# ::mbtools::system_generation::uniform_list --
#
#
proc ::mbtools::system_generation::cylinder::uniform_list_cyl { num range } {

    set step [expr int(floor($range/$num))]

    for {set i 0 } { $i < $num } { incr i } {
	lappend list [expr $i*$step ]
    }
    return $list
}
#
# ::mbtools::system_generation::create_cylinder --
#
# Place a topology onto a cylindrical geometry. This is done by fudging the
# area per molecule until the correct integer number of molecules can
# be placed. It is not advisable to use this for molecules that are of
# very different sizes.
#
# Note: Our cylinder will be periodic in z-direction
#
# Options:
#
# random: place molecules randomly rather than uniformly over the
# cylinder. Useful if you have a very heterogenous system or one that
# just won't work with the uniform method.
#
proc ::mbtools::system_generation::cylinder::create_cylinder { args } {
    mmsg::send [namespace current] "placing lipids on a cylinder"
    # ---- Process Command Line Args --------------------------- #

    set options {
	{center.arg  { 0.0 0.0 0.0 } "location of the sphere relative to the box center" }
	{initarea.arg      1.29    "the starting value for area per mol" }
	{shuffle   "shuffle the topology before placing molecules" }
	{bondl.arg     1.0   "bond length between atoms"  }
    }
    set usage "Usage: create_cylinder \[center:initarea:shuffle:bondl]: "
    array set params [::cmdline::getoptions args $options $usage]

    global ::mbtools::system_generation::moltypeskey
    global ::mbtools::system_generation::boxl
    global ::mbtools::system_generation::topology

    if { $params(shuffle) } {
	set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }

    set center [list 0.0 0.0 0.0 ]
    #Construct the real center
    lset center 0 [expr [lindex $params(center) 0] + [lindex $boxl 0]/(2.0)]
    lset center 1 [expr [lindex $params(center) 1] + [lindex $boxl 1]/(2.0)]
    lset center 2 [expr [lindex $params(center) 2] ]

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

    # Make some aliases for useful numbers
    set del [expr $avmolsize/2.0]
    set pi 3.141592
    set marea $params(initarea)
    set ellz [lindex $boxl 2]		;# length of box in z-direction

    # calculate the radius of our cylinder
    set rad_av [expr $nmols*$marea/(4.0*$pi*$ellz)] ;# average radius of our cylindrical vesicle
    mmsg::send [namespace current] "radius: $rad_av"

    #Based the calculated radius and average molecule size work out
    #how many should go on outer and inner leaflets
    set n_mols_outer [expr 2.0*$pi*$ellz*($rad_av+$del)/$marea]
    set n_mols_inner [expr 2.0*$pi*$ellz*($rad_av-$del)/$marea]
    mmsg::send [namespace current] "n+: $n_mols_outer n-:$n_mols_inner"
    set n_mols_outer [expr int(($n_mols_outer))]
    set n_mols_inner [expr $nmols - $n_mols_outer]
    mmsg::send [namespace current] "n+: $n_mols_outer n-:$n_mols_inner"

    # Determine radius for the outer layer
    set rad [expr $rad_av + $del ]
    mmsg::send [namespace current] "outer radius: $rad"

    # Determine the area per molecule required to place all the lipids
    # in the outer layer.
    set fudgedlipidarea [ find_fudgearea_cyl $rad $boxl $marea $n_mols_outer ]

    #  In some cases no value for the area can be found to give this
    # number of lipids so we then have <exlipids> in excess which
    # should be removed from the sphere.
    set exmols [expr  [find_n_cyl $rad $boxl $fudgedlipidarea ]  - $n_mols_outer]
    if { $exmols > 0 } {
	if { $nmoltypes > 1 } {
	    mmsg::warn [namespace current] "molecule proportions may be incorrect!!"
	}
	set removelist [uniform_list_cyl $exmols [expr $n_mols_outer + $exmols]]
    } else { set removelist [expr $nmols + $exmols + 50]}


    # --- Outer Layer -------- #
    # Now go through and place all mols in the outer layer
    set d [expr sqrt($fudgedlipidarea)]
    set mzeta [expr int(($ellz/$d))]
    set dzeta [expr $ellz/$mzeta]
    set dtheta [expr $fudgedlipidarea/$dzeta]
    set mtheta [expr int((2.0*$pi*$rad/$dtheta))]

    set molnum 0
    set outercount 0
    set removecount 0
    for { set m 0 } { $m < $mzeta } { incr m } {
	set ellz_curr [expr ($m+0.5)*$ellz/($mzeta*1.0)]
	for { set n 0 } { $n < $mtheta } { incr n } {
	    set theta [expr 2.0*$pi*$n/(1.0*$mtheta)]
	    # Find the position of the first tail bead in the molecule
	    lappend tailpos [expr ($rad_av+0.5)*cos($theta) + [lindex $center 0]]
	    lappend tailpos [expr ($rad_av+0.5)*sin($theta) + [lindex $center 1]]
	    lappend tailpos [expr  $ellz_curr               + [lindex $center 2]]

	    # Find the orientation vector for the molecule that is
	    # normal to the sphere
	    lappend orient [expr cos($theta)]
	    lappend orient [expr sin($theta)]
	    lappend orient [expr 1.0        ]

	    if { [lindex $removelist $removecount] != $outercount } {
		# Place the molecule
		set mol [lindex $topology $molnum]
		::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
		incr molnum
		incr outercount
	    } else {
		# Skip placing a lipid at this position
		incr outercount
		incr removecount
	    }
	    unset tailpos
	    unset orient
	}
    }
    mmsg::send [namespace current] "placed [expr $outercount-$removecount] mols in outer layer"
    mmsg::send [namespace current] "$removecount mols were skipped to get this number"

    #------------------ Inner ---------#
    # Determine r1 and r2 for the inner layer
    set rad [expr $rad_av - $del ]
    mmsg::send [namespace current] "inner radius: $rad"

    set fudgedlipidarea [ find_fudgearea_cyl $rad $boxl $marea $n_mols_inner ]
    set exmols [expr  [find_n_cyl $rad $boxl $fudgedlipidarea ]  - $n_mols_inner]
    if { $exmols > 0 } {
	if { $nmoltypes > 1 } {
	    mmsg::warn [namespace current] "molecule proportions may be incorrect!!"
	}
	set removelist [uniform_list_cyl $exmols [expr $n_mols_inner + $exmols]]
    } else { set removelist [expr $nmols + 1]}


    # Now go through and place all mols in the Inner layer
    set d [expr sqrt($fudgedlipidarea)]
    set mzeta [expr int(($ellz/$d))]
    set dzeta [expr $ellz/$mzeta]
    set dtheta [expr $fudgedlipidarea/$dzeta]
    set mtheta [expr int((2.0*$pi*$rad/$dtheta))]

    set innercount 0
    set removecount 0
    for { set m 0 } { $m < $mzeta } { incr m } {
	set ellz_curr [expr ($m+0.5)*$ellz/($mzeta*1.0)]
	for { set n 0 } { $n < $mtheta } { incr n } {
	    set theta [expr 2.0*$pi*$n/(1.0*$mtheta)]
	    lappend tailpos [expr ($rad_av-0.5)*cos($theta) + [lindex $center 0]]
	    lappend tailpos [expr ($rad_av-0.5)*sin($theta) + [lindex $center 1]]
	    lappend tailpos [expr  $ellz_curr               + [lindex $center 2]]

	    # Find the normal vector pointing inwards
	    lappend orient [expr -cos($theta)]
	    lappend orient [expr -sin($theta)]
	    lappend orient [expr -1.0]

	    if { [lindex $removelist $removecount] != $innercount } {
		# Place the molecule
		set mol [lindex $topology $molnum]
		::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
		incr molnum
		incr innercount
	    } else {
		incr innercount
		incr removecount
	    }
	    unset tailpos
	    unset orient

	}
    }


     # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }

    mmsg::send [namespace current] "placed [expr $innercount-$removecount] mols in inner layer"
    mmsg::send [namespace current] "$removecount random mols were skipped to get this number"
    mmsg::send [namespace current] "total mols in vesicle: $molnum"
    mmsg::send [namespace current] "uniform cylindrical vesicle created"
    flush stdout

    return
}
