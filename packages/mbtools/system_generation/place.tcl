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
#
# Routines for placing molecules
#



namespace eval mbtools::system_generation {}

namespace eval ::mbtools::system_generation::placemol {
    #variable dists 
}



# ::mbtools::system_generation::placemol-- 
#
# general routine for placing molecules
#
#
proc ::mbtools::system_generation::placemol { mol pos args } { 
    set options {
	{bondl.arg     1.0   "bond length between atoms"  }
	{orient.arg  { 0 0 1 } "orientation vector for the mol " }
	{orient2.arg  { 1 0 0 } "a second orientation vector for the mol " }

    }
    set usage "Usage: create_bilayer topo boxl \[bondl:orient]"
    array set params [::cmdline::getoptions args $options $usage]
    

    # Retrieve the molecule information for this molecule type	
    set typekey [matchtype [lindex $mol 0]]

    # Place the molecule depending on type
    switch [lindex $typekey 1] {
	"lipid" {
	    place_lipid $mol $params(orient) $pos -bondl $params(bondl) 
	}
	"hollowsphere" {
	    place_hollowsphere $mol $pos $params(orient)
	}
	"spanlipid" {
	    place_lipid $mol $params(orient) $pos -bondl $params(bondl) -midpos 
	}
	"sphericalcap" {
	    place_sphericalcap $mol  $pos $params(orient)
	}
	"cylindricalsegment" {
	    place_cylindricalsegment $mol  $pos $params(orient) $params(orient2)
	}
	"hexagon" {
	    place_hexagon $mol  $pos $params(orient)
	}

	"protein" {
	    place_protein $mol $pos $params(orient)
	}
	"sphericalconstraint" {
	    place_sphericalconstraint $mol $pos
	}
	"default" {
	    ::mmsg::err [namespace current] "couldn't place molecule of type [lindex $typekey 1], possibilities are: \n lipid \n hollowsphere \n sphericalconstraint"
	}
    }

    return
}

# ::mbtools::system_generation::matchtype-- 
#
# Search the molecule type list for a type number and return the type specification
#
#
proc ::mbtools::system_generation::matchtype { mol } {

    variable moltypeskey

    foreach key $moltypeskey {
	if { [lindex $key 0] == [lindex $mol 0] } {
	    return $key
	}
    }
    mmsg::err [namespace current] "could not find a matching key to moltype [lindex $mol 0]"
}




# ::mbtools::system_generation::place_bonds--
#
# place bonds between particles.  adjust lengths of bonds appropriately
# usefull routine for creating stiff structures
#
# mol : list of atoms in the molecule
# cutoff : cutoff distance for placing bonds
# params : lj params for repulsions between the molecules (only between first half of atoms)
# dists  : a list of the lengths of the bonds
# startingbondno : index number of first bond in molecule - they then go sequentially
proc ::mbtools::system_generation::place_bonds {mol cutoff params dists startingbondno minpartnum laynum} {
    set harm_k [lindex $params 0]
    set lj_eps [lindex $params 1]
    set lj_sigma [lindex $params 2]
    set lj_cutoff [lindex $params 3]
    # Figure out which beads are neigbours (2*$midist) is used to decide who should be bonded to who)
    set totnum [expr ([llength $mol ] -1) ]
    for { set i 1 } { $i <= $totnum } { incr i } {
	set partnum [lindex $mol  $i ]
	lappend nbrs [analyze nbhood $partnum $cutoff]
    }

    puts "dists at beginning are $dists"
    
    # Bond the neighbours
    # We require a unique type of bond for each length of bond required.
    # Thus we work out how many types of bonds we need, of what lengths, and between which particles they should be applied.
    # The bond length takes into account the lj repulsions between the particles such that the resting ground state length of the bond
    # corresponds to what we geometrically require for our spherical cap.
    set count [llength $dists]
    for { set i 1 } { $i <= $totnum } { incr i } {
	set partnum [lindex $mol  $i ]
	# get a list of the neighbours for this particle
	set nblist [lindex $nbrs [expr $i - 1]]
	# find out what the distance is to each of its neighbours
	foreach atom $nblist {	    
	    if { $atom > $partnum } {
		set part1pos [part $partnum print pos]
		set part2pos [part $atom print pos]
		set boxl [setmd box_l]
		set xdist [expr [lindex $part1pos 0]-[lindex $part2pos 0]]
		set ydist [expr [lindex $part1pos 1]-[lindex $part2pos 1]]
		set zdist [expr [lindex $part1pos 2]-[lindex $part2pos 2]]
		set box_l [setmd box_l]
		while {$xdist > [expr [lindex $box_l 0]/2.0]} {
		    set xdist [expr $xdist - [lindex $box_l 0]]
		}
		while {$xdist < [expr -[lindex $box_l 0]/2.0]} {
		    set xdist [expr $xdist + [lindex $box_l 0]]
		}
		while {$ydist > [expr [lindex $box_l 1]/2.0]} {
		    set ydist [expr $ydist - [lindex $box_l 1]]
		}
		while {$ydist < [expr -[lindex $box_l 1]/2.0]} {
		    set ydist [expr $ydist + [lindex $box_l 1]]
		}
		while {$zdist > [expr [lindex $box_l 2]/2.0]} {
		    set zdist [expr $zdist - [lindex $box_l 2]]
		}
		while {$zdist < [expr -[lindex $box_l 2]/2.0]} {
		    set zdist [expr $zdist + [lindex $box_l 2]]
		}
		set dist [expr pow($xdist,2) +  pow($ydist,2) + pow($zdist,2)]
		set dist [expr pow($dist,0.5)]
		set exist 0
		# have we already created a bond with this length
		#puts "dists is $dists"
		for { set j 0} {$j < $count} {incr j} {
		    set val [lindex $dists $j]
		    if {[expr pow($dist - $val,2)] < 0.00000001} {
			set exist 1
			set bondno [expr $j + $startingbondno]
			# we have already created one, so just apply another instance of it
			#puts "setting bond between $partnum and $atom of type $bondno"
			#puts "part $partnum bond $bondno $atom with dist $dist"
			part $partnum bond $bondno $atom
			set j $count
		    }
		}
		# if we haven't already created one of this length then
		if {$exist == 0} {
		    lappend dists $dist
		    set bondno [expr $count + $startingbondno]
		    incr count
		    set actualdist $dist
		    # do they two atoms have a nonbonded forces between them
		    # if they do then compensate for it in the bond length
		    if {($partnum <= $minpartnum + $laynum - 1) && ($atom <= $minpartnum + $laynum -1) && ($dist < $lj_cutoff)} {
			set lj_force [expr - 24.0 * $lj_eps / $dist * (2.0 * pow($lj_sigma/$dist,12)-pow($lj_sigma/$dist,6))]
			if {$lj_force > 0} {
			    puts "****************OH DEAR --------- error is setting up sphere "
			}
			set harm_offset [expr pow(-$lj_force/$harm_k,0.5)]
			set dist [expr $dist - $harm_offset]
		    }
		    # create new bond type and apply it to pair
		    inter $bondno harmonic $harm_k $dist
		    puts "creating new bond number $bondno with length $actualdist"
		    #puts "part $partnum bond $bondno $atom with dist $actualdist"
		    part $partnum bond $bondno $atom 
		}
	    }
	}
    }
    return $dists
}









# ::mbtools::system_generation::place_lipid-- 
#
# Place a lipid with all atoms in a line.
#
# Arguments:
#
# orient: A vector along which to place the lipids
# mol: The molecule to be placed
# pos: The position of the first tail bead
#
# Note that both linkbond and bendbond will have previously been set by set_bonded_interactions
#
# Assumes head bead is given first
#
#
proc ::mbtools::system_generation::place_lipid { mol orient pos args } {
    set options {
	{bondl.arg     1.0   bond length between atoms  }
	{midpos  "the value of pos corresponds to the middle bead not the tail"}
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:uniform]"
    array set params [::cmdline::getoptions args $options $usage]

    # Ensure that the orientation vector is normalized
    set orient [::mbtools::utils::normalize $orient]


    # Max a few  aliases
    set rx [lindex $pos 0]
    set ry [lindex $pos 1]
    set rz [lindex $pos 2]
    set nx [lindex $orient 0]
    set ny [lindex $orient 1]
    set nz [lindex $orient 2]

    set moltype [lindex $mol 0]
    set typeinfo [matchtype $moltype]
 
    set nbeads [expr [llength $mol] -1]

    # Determine the lipid length and enforce the midpos condition if it
    # exists
    set lipidlen [expr $nbeads*$params(bondl)]
    set halflen [expr $lipidlen/2.0]
    if { $params(midpos) } {
	set rx [expr $rx - $halflen*$nx]
	set ry [expr $ry - $halflen*$ny]
	set rz [expr $rz - $halflen*$nz]
    }
    # -- Extract the bonding information -- #

    set bond_params [lindex $typeinfo 3 ]


    for { set b 0 } { $b < $nbeads } {incr b } {
	set partnum [lindex $mol [expr $b + 1]]

	# Note that because the head bead is given first but we start
	# placement     variable av_molcom_ifrom the tail we need to use the following
	# formulas for position
	set posx [expr $rx+($nbeads - $b -1.0)*$params(bondl)*$nx]
	set posy [expr $ry+($nbeads - $b -1.0)*$params(bondl)*$ny]
	set posz [expr $rz+($nbeads - $b -1.0)*$params(bondl)*$nz]
	
	# -- Determine the particle type  ----
	set ptype [lindex $typeinfo 2 $b]


	part $partnum pos $posx $posy $posz type $ptype 

	#As soon as the first particle has been placed we can start
	#connecting bonds
	if { $b > 0 } {
	    set a1 [lindex $mol $b]
	    set a2 $partnum
	    set bt [lindex $bond_params 0 ]
	    
	    part $a1 bond $bt $a2
	}
	
	# And placing Pseudo bending potentials
	if { $b > 1 } {
	    set a1 [lindex $mol [expr $b -1] ]
	    set a2 $partnum
	    set bt [lindex $bond_params 1 ]
	    part $a1 bond $bt $a2
	}	
    }


}


# ::mbtools::system_generation::place_sphericalconstraint-- 
#
# Place a spherical constraint
#
# Arguments:
#
# mol: Information on the particle type for the constraint
# pos: The center of the sphere
#
proc ::mbtools::system_generation::place_sphericalconstraint { mol pos args } {
    set options {
    }
    set usage "Usage: place_sphericalconstraint topo boxl \[mol:pos:]"
    array set params [::cmdline::getoptions args $options $usage]


    set ptype [lindex $mol 0]
    set typekey [matchtype $ptype]

    set radius [lindex $typekey  3 0]
    set direction [lindex $typekey  3 1]

    constraint sphere center [lindex $pos 0] [lindex $pos 1] [lindex $pos 2] radius $radius direction $direction type $ptype

    return

}



#
# Arguments:
#
# mol: particle types and the molecule type id
# pos: The beads' positions
#
proc ::mbtools::system_generation::place_protein { mol pos orient args } {

    variable middlebead

    set moltype [lindex $mol 0]
    set typeinfo [matchtype $moltype]

    set atomtypes [lindex $typeinfo 2]
    set bonds [lindex $typeinfo 3]

  # Ensure that the orientation vector is normalized
    set orient [::mbtools::utils::normalize $orient]

    #set bond_params [lindex $moltypeskey $moltype 3]

    #nbeads = number of beads in a protein

    set nbeads [expr [llength $mol] -1]

    set rx [lindex $pos 0]
    set ry [lindex $pos 1]
    set rz [lindex $pos 2]

    set nx [lindex $orient 0]
    set ny [lindex $orient 1]
    set nz [lindex $orient 2]
   
    # Create a protein with 2 hexagonal layers and has 12 layer from bottom to up#
    # Placing the beads # # $a is the height of the protein #

    for {set a 0 } { $a < 12 } {incr a} {

    set c [expr (19*$a)]
    set d [expr (19+19*$a)]
    set e [expr (7+19*$a)]
    set f [expr (1+19*$a)]

	 for { set b $c } { $b < $d } {incr b} {

	     if { $b == $c } {

				     
		    set partnum [lindex $mol [expr $b +1]] 
		    set ptype [lindex $typeinfo 2 $b]

		    set posx [expr $rx + ($a*1.0)*$nx]
		    set posy [expr $ry + ($a*1.0)*$ny]
		    set posz [expr $rz + ($a*1.0)*$nz]

  
		    part $partnum pos $posx $posy $posz type $ptype
		
		} else { if { $b < $e } {
		     
		    set partnum1 [lindex $mol [expr $b + 1]]
		    set ptype [lindex $typeinfo 2 $b] 
	   
		  
	  	    set posx [expr $rx+ (0.9+0.1*$a)*cos($b*(2*3.14)/(6))] 
		    set posy [expr $ry+ (0.9+0.1*$a)*sin($b*(2*3.14)/(6))]

		    set posz [expr $rz + ($a*1.0)*$nz]
		  
	  	   
		    part $partnum1 pos $posx $posy $posz type $ptype
		
		} else {
		   
		    set partnum1 [lindex $mol [expr $b + 1]]
		    set ptype [lindex $typeinfo 2 $b] 
	   
		  
	  	    set posx [expr $rx+ 2*(0.9+(0.1*$a))*cos($b*(2*3.14)/(12))] 
		    set posy [expr $ry+ 2*(0.9+(0.1*$a))*sin($b*(2*3.14)/(12))]

		    set posz [expr $rz + ($a*1.0)*$nz]
		    
	  	   
		    part $partnum1 pos $posx $posy $posz type $ptype
		} 
	
		}
 
	 }
  
}


    #Connecting beads with bonds#

    for {set a 0 } { $a < 12 } {incr a} {


	set d [expr (19+19*$a)]

	set f [expr (1+19*$a)]

	for { set i $f } { $i <= $d } { incr i } {
	
	set partnum1 [lindex $mol  $i ]
	lappend nbrs [analyze nbhood $partnum1 [expr 0.6 + (0.9+(0.1*$a))]]
	
	set ba [lindex $bonds $a]

	set nblist [lindex $nbrs [expr $i - 1]]

	foreach atom $nblist {	    
	    if { $atom > $partnum1 } {
#		Only place the bond if we haven't already done so.

		part $partnum1 bond $ba $atom
	
	    }
	}
    }
}


#Get the midlle bead id#

 
    for { set b 95 } { $b < 96 } {incr b} {

   
	lappend middlebead [lindex $mol [expr $b+1]]


    } 
  
  

return
}




# ::mbtools::system_generation::place_sphericalcap-- 
#
# Construct a spherical cap
#
#

# Note that this routine uses the icosahedral codes from R. H. Hardin,
# N. J. A. Sloane and W. D. Smith .  In order to get this to work you
# first need to make sure that you have their program creconstruct in
# your path and that their script icover.sh runs fine from the command
# line in the directory where you execute espresso.  You then need to
# ensure that the number of atoms you used corresponds to one which
# they have actually tabulated.

#
# Arguments:
#
# mol: particle types and the molecule type id
# pos: The center of the sphere
#
proc ::mbtools::system_generation::place_sphericalcap { mol pos orient } {
    variable icovermagicnums
    variable dists

    # check to see if dists has already been worked out.  i.e. has a hollowsphere already been placed?
    set err [catch {set dists $dists}]
    if {$err == 1} {
	set dists ""
    } 

    puts "orientation is $orient"

    # what molecule type (number) is this
    set moltype [lindex $mol 0]
    # info on molecules of this type
    set typeinfo [matchtype $moltype]
    # A number of bonds may be defined in the process of creating this molecule.
    # this variable defines the index number of the first such bond.  Other follow sequentially.
    set startingbondno [lindex $typeinfo 3]
    # the radii of the inner and outer shells of the sphere
    set radii [lindex $typeinfo 4]
    set rout [lindex $radii 0]
    set rin [lindex $radii 1]
    # the number of atoms it would take to cover the entire sphere if we were generating the entire sphere
    set natomscov [lindex $typeinfo 5]
    # the lj parameters for the nb interactions between atoms in the molecule (must all be the same)
    set int_params [lindex $typeinfo 6]
    set harm_k [lindex $int_params 0]
    set lj_eps [lindex $int_params 1]
    set lj_sigma [lindex $int_params 2]
    set lj_cutoff [lindex $int_params 3]
    # the total number of atoms in our spherical cap
    set totnum [expr ([llength $mol ] -1) ]
    # the total number of atoms in each layer
    set laynum [expr $totnum / 2]

    #---------------------------------------------------------------------------------------------------------------
    # use the icover routine to get a list of coordinates for particles on a sphere
    set imagic [lsearch -integer -exact $icovermagicnums $natomscov]
    if { $imagic == -1 } {
	foreach val $icovermagicnums {
	    if { $val > $natomscov } {
		set isuggest [expr $imagic ]
		break;
	    }
	    incr imagic
	}
	mmsg::err [namespace current] "can't construct hollowsphere because $natomscov is not an icover magic number try [lindex $icovermagicnums $isuggest] or [lindex $icovermagicnums [expr $isuggest + 1]]"
    }
    if { [catch { set cov [exec icover.sh $natomscov  ] } ] } {
	::mmsg::err [namespace current] "couldn't construct hollowsphere because errors occured when trying to run icover with $natomscov atoms.  icover extracts the icosahedral codes which are copyright R. H. Hardin, N. J. A. Sloane and W. D. Smith, 1994, 2000. so you should obtain them yourself from http://www.research.att.com/~njas/"
    } else {


	set ncov [expr int([llength $cov]/3.0)]

	if { $natomscov != $ncov } {
	    mmsg::err [namespace current] "icover.sh returned $ncov atoms but our hollowsphere has $natomscov"
	}

	# Now sort all the data in cov into a list of points
	for { set i 0 } { $i < $ncov } { incr i } {
	    lappend tmp [lindex $cov [expr 3*$i]]
	    lappend tmp [lindex $cov [expr 3*$i + 1]]
	    lappend tmp [lindex $cov [expr 3*$i + 2]]
	    lappend coords $tmp
	    unset tmp
	}
    }
    # -------------------------------------------------------------------------------------------------------------
    # neighbours is number of atoms withing a certain distance "cutoff" of each atom
    # we use neighbours to orientate the spherical cap (i.e. we want a centre of symmetry to correspond to the centre of our spherical cap)
    for {set i 0} {$i < $ncov} {incr i} {
	lappend neighbours 0
    }
    set cutoff [expr 1.5*sqrt(12/($ncov*1.0))]
    
    # caculate neighbours for each atom
    set mdist -1
    for {set i 0} {$i < $ncov} {incr i} {
	for {set j [expr $i + 1]} {$j < $ncov} {incr j} {
	    set distance [::mbtools::utils::distance [lindex $coords $i] [lindex $coords $j]]
	    if {$distance  < $cutoff} {
		set neighbours [lreplace $neighbours $i $i [expr [lindex $neighbours $i] + 1]]
		set neighbours [lreplace $neighbours $j $j [expr [lindex $neighbours $j] + 1]]
		if {($distance < $mdist) || ($mdist == -1)} {
		    set mdist $distance
		}
	    }
	}
    }

    # numneigh_freq givens the number of atoms which have a certain number of neighbours
    set jcount 0
    for {set i 0} {$i < $ncov} {incr i} {
	set alreadythere 0
	for {set j 0} {$j < $jcount} {incr j} {
	    if {[lindex $numneigh $j] == [lindex $neighbours $i]} {
		set alreadythere 1
		set numneigh_freq [lreplace $numneigh_freq $j $j [expr 1+[lindex $numneigh_freq $j] ]]
		set j $jcount
	    }
	}
	if {$alreadythere == 0} {
	    lappend numneigh [lindex $neighbours $i]
	    lappend numneigh_freq 1
	    incr jcount
	}
    }
    set minfreq -1
    # we then look to see what is the rarest number of neighbours to have
    for {set j 0} {$j < $jcount} {incr j} {
	if {($minfreq == -1) || ($minfreq > [lindex $numneigh_freq $j])} {
	    set minfreq [lindex $numneigh_freq $j]
	    set numneigh_rarest [lindex $numneigh $j]
	}
    }
    # we select an atom from this rare group.  This will be the centre of our spherical cap.
    for {set i 0} {$i < $ncov} {incr i} {
	if {[lindex $neighbours $i] == $numneigh_rarest} {
	    set rare_atom $i
	    set i $ncov
	}
    }

    set coord_centre [lindex $coords $rare_atom]


    # we know select the $laynum many atoms closest to the rare atom we have selected as the centre.
    # previously we had coordinates for an entire sphere. 
    set coordslength $ncov
    for {set i 0} {$i < $laynum} {incr i} {
	set minind -1
	for {set j 0} {$j < $coordslength} {incr j} {
	    set val [::mbtools::utils::distance $coord_centre [lindex $coords $j]]
	    if {($minind == -1) || ($val < $minval)} {
		set minind $j
		set minval $val    
		
	    }
	}
	lappend newcoords [lindex $coords $minind]
	set coords [lreplace $coords $minind $minind]
	set coordslength [expr $coordslength - 1]
    }

    set coords $newcoords

    for {set i 0} {$i < $laynum} {incr i} {
	if { [lindex $coords $i] == $coord_centre} {
	    set rare_atom $i
	}
    }
    

    # we rotate the coordinates of the particles in the sphere such that the chosen atom lies at the bottom of the sphere (minimum z value).
    set coord_centre [::mbtools::utils::normalize $coord_centre]
    set orient [::mbtools::utils::normalize $orient]
    set desired_centre [::mbtools::utils::scalevec $orient -1]
    set rotation_angle [expr acos([::mbtools::utils::dot_product $coord_centre $desired_centre])]
    set rotation_axis [::mbtools::utils::cross_product $coord_centre $desired_centre]
    set rotation_axis [::mbtools::utils::normalize $rotation_axis]
    set rotation_matrix [::mbtools::utils::rotation_matrix $rotation_axis $rotation_angle]

    set tmp ""
    foreach coord $coords {
	lappend tmp [::mbtools::utils::matrix_vec_multiply $rotation_matrix $coord]
    }
    set coords $tmp

    set minpartnum -1

    # Place the beads in preliminary positions on the outer surface of the spherical cap
    for { set i 0 } { $i < $laynum } { incr i } {
	set tmp [lindex $coords $i]
	set partnum [lindex $mol  [expr $i +1]]
	set parttype [lindex $typeinfo 2 $i]
	set xpos [expr $rout * [lindex $tmp 0] + [lindex $pos 0]]
	set ypos [expr $rout * [lindex $tmp 1] + [lindex $pos 1]]
	set zpos [expr $rout * [lindex $tmp 2] + [lindex $pos 2]]
	part $partnum pos $xpos $ypos $zpos type $parttype
	if {($partnum < $minpartnum) || ($minpartnum == -1)} {
	    set minpartnum $partnum
	}
    }    

    
    set mdist [expr $mdist * $rout]

    # Place the beads in preliminary positions on the inner surface of the spherical cap
    for { set i 0 } { $i < $laynum } { incr i } {
	set tmp [lindex $coords $i]
	set partnum [lindex $mol  [expr $i +1 + $laynum]]
	set parttype [lindex $typeinfo 2 [expr $i + $laynum]]
	set xpos [expr $rin * [lindex $tmp 0] + [lindex $pos 0]]
	set ypos [expr $rin * [lindex $tmp 1] + [lindex $pos 1]]
	set zpos [expr $rin * [lindex $tmp 2] + [lindex $pos 2]]
	part $partnum pos $xpos $ypos $zpos type $parttype
	if {($partnum < $minpartnum) || ($minpartnum == -1)} {
	    set minpartnum $partnum
	}
    }    

    # Based on a desired value of mdist equal to 1.0 find the required radius
    mmsg::send [namespace current] "creating hollow sphere with radius $rout, $rin and $natomscov beads on both shells"

    set cutoff [expr $mdist * 2.0]

    set dists [::mbtools::system_generation::place_bonds $mol $cutoff $int_params $dists $startingbondno $minpartnum $laynum]

    return

}


proc ::mbtools::system_generation::place_cylindricalsegment { mol pos orient orient2} {
    variable dists
    
    # check to see if dists has already been worked out.  i.e. has a hollowsphere already been placed?
    set err [catch {set dists $dists}]
    if {$err == 1} {
	set dists ""
    } 

    # what molecule type (number) is this
    set moltype [lindex $mol 0]
    # info on molecules of this type
    set typeinfo [matchtype $moltype]
    # A number of bonds may be defined in the process of creating this molecule.
    # this variable defines the index number of the first such bond.  Other follow sequentially.
    set startingbondno [lindex $typeinfo 3]
    # the radii of the inner and outer shells of the sphere
    set radii [lindex $typeinfo 4]
    set rout [lindex $radii 0]
    set rin [lindex $radii 1]
    # the lj parameters for the nb interactions between atoms in the molecule (must all be the same)
    set int_params [lindex $typeinfo 5]
    set harm_k [lindex $int_params 0]
    set lj_eps [lindex $int_params 1]
    set lj_sigma [lindex $int_params 2]
    set lj_cutoff [lindex $int_params 3]
    # the total number of atoms in our spherical cap
    set totnum [expr ([llength $mol ] -1) ]
    # the total number of atoms in each layer
    set laynum [expr $totnum / 2]
    # the surface of the cylinder is covered in a hexagonal lattice of points
    set shape [lindex $typeinfo 6]
    # number of layers in direction1
    set width1 [lindex $shape 0]    
    # number of layers in direction2
    set width2 [lindex $shape 1]
    # bond length
    set bl [lindex $shape 2]

    # since the beads are placed on a hexagon lattice the two sets of edges will be different from each other
    # along two edges the beads will all be in a straight line
    # along the other the beads will go in a zigzag
    # for this reason it is necessary to specify in which direction the cylinder should be curved and
    # whether the zigzag edges should be symmetric or antisymmetric.

    # the direction in which the cylinder should be curved
    # the options are 1,0,-1 (0 corresponds to no curvature)
    set curvedir [lindex $shape 3]

    # whether the edges shuld be symmetric or antisymmetric
    # set edgesim to 0 unless you want cylinder to span the box in which case set it to 1
    set edgesym [lindex $shape 4]

    # calculate the total number of bead required to make the cylinder based on shape parameters given
    if {$edgesym == 1} {
	set totnum_implied [expr $width2*$width1*2.0]
    } else {
	set totnum_implied [expr (($width2 + 0.5) * $width1 -0.5)*2.0 ]
	set edgesym 0
    }

    # compare the beads required to the list of bead given.
    if {$totnum_implied != $totnum} {
	puts "the molcule parameters are not consistent. according to the two widths given ($width1 , $width2) there should be $totnum_implied atoms but there are $totnum in the specifications."
	return 1
    }

    if {$curvedir == 0} {
	puts "curvedir is 0 so molecule will not be curved"
    }

    #calculate the rotation matrix necessary to rotate the cylinder into the correct position
    set initial_orient {0 0 1}
    set initial_orient2 {1 0 0}
    set orient [::mbtools::utils::normalize $orient]
    set rotation_angle [expr acos([::mbtools::utils::dot_product $initial_orient $orient])]
    set rotation_axis [::mbtools::utils::cross_product $initial_orient $orient]
    if {[::mbtools::utils::distance $rotation_axis {0 0 0}] == 0} {
	# if initial_orient and orient are parellel then we  can use any vector perpendicular to this line as the rotation axis
	set rotation_axis [::mbtools::utils::perp_vec $orient]
    }
    set rotation_axis [::mbtools::utils::normalize $rotation_axis]
    set rotation_matrix [::mbtools::utils::rotation_matrix $rotation_axis $rotation_angle]

    set initial_orient2 [::mbtools::utils::matrix_vec_multiply $rotation_matrix $initial_orient2]
    set orient2 [::mbtools::utils::scalevec $orient2 [expr 1 - [::mbtools::utils::dot_product $orient2 $orient]]]
    set orient2 [::mbtools::utils::normalize $orient2]
    set rotation_angle [expr acos([::mbtools::utils::dot_product $initial_orient2 $orient2])]
    set rotation_axis [::mbtools::utils::cross_product $initial_orient2 $orient2]
    if {[::mbtools::utils::distance $rotation_axis {0 0 0}] == 0} {
	# if initial_orient and orient are parellel then we  can use any vector perpendicular to this line as the rotation axis
	set rotation_axis [::mbtools::utils::perp_vec $orient2]
    }
    set rotation_axis [::mbtools::utils::normalize $rotation_axis]
    set rotation_matrix2 [::mbtools::utils::rotation_matrix $rotation_axis $rotation_angle]

    set rotation_matrix [::mbtools::utils::matrix_multiply $rotation_matrix2 $rotation_matrix]

    # loop over the two widths and place the particles
    set atom 0
    for {set i 0} {$i < $width1} {incr i} {
	if {$i%2 == 0} {
	    set thick $width2
	    # the variable 'adjust' is used to make the edges of the cylinder symmetric or antisymmetric
	    set adjust 0
	} else {
	    set thick [expr $width2+1]
	    set adjust [expr $edgesym]
	}
	for {set j 0} {$j < [expr $thick - $adjust]} {incr j} {
	    #initially set x and y to positions on a rectangle
	    set x [expr (-($width1-1.0)/2.0 + $i)*$bl*pow(3,0.5)/2.0]
	    set y [expr (-($thick-1.0)/2.0 + $j)*$bl]

	    #place the particles in the outer layer
	    if {$curvedir == 0} {
		set z 0
	    } else  {
		# then reset the x and y coordinates and set the z coordinate so that the rectangle curves upwards
		if {$curvedir > 0 } {
		    set arcl $y
		} else {
		    set arcl $x
		}
		set angle [expr $arcl/$rout]
		set xy [expr $rout * sin($angle)]	    
		set z [expr $rout * (1-cos($angle))]
		if {$curvedir > 0 } {
		    set y $xy
		} else {
		    set x $xy
		}
	    }
	    set rotated_coord [::mbtools::utils::matrix_vec_multiply $rotation_matrix "$x $y $z"]
	    set x [lindex $rotated_coord 0]
	    set y [lindex $rotated_coord 1]
	    set z [lindex $rotated_coord 2]
	    # set the rectangle to have the correct base point
	    set x [expr $x + [lindex $pos 0]]
	    set y [expr $y + [lindex $pos 1]]
	    set z [expr $z + [lindex $pos 2]]
	    set partnum [lindex $mol  [expr $atom +1]]
	    set parttype [lindex $typeinfo 2 $atom]
	    part $partnum pos $x $y $z type $parttype

	    #place the particles in the inner layer
	    if {$curvedir == 0} {
		set z 0
	    } else {
		set xy [expr $rin * sin($angle)]	    
		set z [expr $rout - $rin*cos($angle)]
		if {$curvedir > 0 } {
		    set y $xy
		    set x [expr (-($width1-1.0)/2.0 + $i)*$bl*pow(3,0.5)/2.0]
		} else {
		    set x $xy
		    set y [expr (-($thick-1.0)/2.0 + $j)*$bl]
		}
	    }
	    set rotated_coord [::mbtools::utils::matrix_vec_multiply $rotation_matrix "$x $y $z"]
	    set x [lindex $rotated_coord 0]
	    set y [lindex $rotated_coord 1]
	    set z [lindex $rotated_coord 2]
	    set x [expr $x + [lindex $pos 0]]
	    set y [expr $y + [lindex $pos 1]]
	    set z [expr $z + [lindex $pos 2]]
	    set partnum [lindex $mol  [expr $atom +$laynum +1]]
	    set parttype [lindex $typeinfo 2 [expr $atom+$laynum]]
	    part $partnum pos $x $y $z type $parttype

	    incr atom
	}
    }
    set cutoff [expr $bl * 2.0]
    set minpartnum [lindex $mol 1]

    set dists [::mbtools::system_generation::place_bonds $mol $cutoff $int_params $dists $startingbondno $minpartnum $laynum]
    puts "dists are $dists"
   
}



proc ::mbtools::system_generation::place_hexagon { mol pos orient } {
    variable dists
    
    # check to see if dists has already been worked out.  i.e. has a hollowsphere already been placed?
    set err [catch {set dists $dists}]
    if {$err == 1} {
	set dists ""
    } 

    # what molecule type (number) is this
    set moltype [lindex $mol 0]
    # info on molecules of this type
    set typeinfo [matchtype $moltype]
    # A number of bonds may be defined in the process of creating this molecule.
    # this variable defines the index number of the first such bond.  Other follow sequentially.
    set startingbondno [lindex $typeinfo 3]
    # the lj parameters for the nb interactions between atoms in the molecule (must all be the same)
    set int_params [lindex $typeinfo 4]
    set harm_k [lindex $int_params 0]
    set lj_eps [lindex $int_params 1]
    set lj_sigma [lindex $int_params 2]
    set lj_cutoff [lindex $int_params 3]
    # the total number of atoms in our spherical cap
    set totnum [expr ([llength $mol ] -1) ]
    # the total number of atoms in each layer
    set laynum [expr $totnum / 2]
    # the surface of the cylinder is covered in a hexagonal lattice of points
    # number of layers in direction1
    set bl [lindex $typeinfo 5]

    set length 0
    set x [lindex $pos 0]
    set y [lindex $pos 1]
    set z [lindex $pos 2]
    set Dy [expr $bl * sqrt(3.0) / 2.0]
    set Dx [expr $bl/2.0]
    set atom 0
    if {$atom < $laynum} {
	set partnum [lindex $mol  [expr $atom +1]]
	set parttype [lindex $typeinfo 2 $atom]
	part $partnum pos $x $y $z type $parttype
	set partnum [lindex $mol  [expr $atom +$laynum +1]]
	set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
	part $partnum pos $x $y [expr $z+$bl] type $parttype
    }
    incr atom
    while {$atom < $laynum} {
	incr length
	for {set i 0} {$i < $length} {incr i} {
	    set x [expr $x + $Dx]
	    set y [expr $y + $Dy]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
	for {set i 0} {$i < [expr $length - 1.0]} {incr i} {
	    set x [expr $x - $Dx]
	    set y [expr $y + $Dy]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
	for {set i 0} {$i < [expr $length]} {incr i} {
	    set x [expr $x - $bl]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
	for {set i 0} {$i < [expr $length]} {incr i} {
	    set x [expr $x - $Dx]
	    set y [expr $y - $Dy]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
	for {set i 0} {$i < [expr $length]} {incr i} {
	    set x [expr $x + $Dx]
	    set y [expr $y - $Dy]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
	for {set i 0} {$i < [expr $length]} {incr i} {
	    set x [expr $x + $bl]
	    if {$atom < $laynum} {
		set partnum [lindex $mol  [expr $atom +1]]
		set parttype [lindex $typeinfo 2 $atom]
		part $partnum pos $x $y $z type $parttype
		set partnum [lindex $mol  [expr $atom +$laynum +1]]
		set parttype [lindex $typeinfo 2 [expr $atom + $laynum]]
		part $partnum pos $x $y [expr $z+$bl] type $parttype
	    }
	    incr atom
	}
    }

    set cutoff [expr $bl * 2.0]
    set minpartnum [lindex $mol 1]

    set dists [::mbtools::system_generation::place_bonds $mol $cutoff $int_params $dists $startingbondno $minpartnum $laynum]
   
}


# ::mbtools::system_generation::place_hollowsphere-- 
#
# Construct a large hollow sphere from small beads
#
#

# Note that this routine uses the icosahedral codes from R. H. Hardin,
# N. J. A. Sloane and W. D. Smith .  In order to get this to work you
# first need to make sure that you have their program creconstruct in
# your path and that their script icover.sh runs fine from the command
# line in the directory where you execute espresso.  You then need to
# ensure that the number of atoms you used corresponds to one which
# they have actually tabulated.

# Another noteworthy point is that you might need to fiddle a bit with
# the cutoff for determining the nearest neighbours.  If you get bonds
# broken messages then quite likely you need to decrease it.  If you
# get a mesh with holes then you need to increase it.

# Assume the last atom type are the filler beads and place the
# first atom types proportionately and in separate phases on the
# sphere

#
# Arguments:
#
# mol: particle types and the molecule type id
# pos: The center of the sphere
#
proc ::mbtools::system_generation::place_hollowsphere { mol pos orient } {
    variable icovermagicnums
    set moltype [lindex $mol 0]
    set typeinfo [matchtype $moltype]

    set atomtypes [lindex $typeinfo 2]
    set bonds [lindex $typeinfo 3]
    set nfill [lindex $typeinfo 4]

    set natomscov [expr [llength $mol ] -1 -$nfill ]


    set imagic [lsearch -integer -exact $icovermagicnums $natomscov]
    if { $imagic == -1 } {
	foreach val $icovermagicnums {
	    if { $val > $natomscov } {
		set isuggest [expr $imagic ]
		break;
	    }
	    incr imagic
	}
	mmsg::err [namespace current] "can't construct hollowsphere because $natomscov is not an icover magic number try [lindex $icovermagicnums $isuggest] or [lindex $icovermagicnums [expr $isuggest + 1]]"
    }

    if { [catch { set cov [exec icover.sh $natomscov  ] } ] } {
	::mmsg::err [namespace current] "couldn't construct hollowsphere because errors occured when trying to run icover with $natomscov atoms.  icover extracts the icosahedral codes which are copyright R. H. Hardin, N. J. A. Sloane and W. D. Smith, 1994, 2000. so you should obtain them yourself from http://www.research.att.com/~njas/"
    } else {


	set ncov [expr int([llength $cov]/3.0)]

	if { $natomscov != $ncov } {
	    mmsg::err [namespace current] "icover.sh returned $ncov atoms but our hollowsphere has $natomscov"
	}

	# Now sort all the data in cov into a list of points
	for { set i 0 } { $i < $ncov } { incr i } {
	    lappend tmp [lindex $cov [expr 3*$i]]
	    lappend tmp [lindex $cov [expr 3*$i + 1]]
	    lappend tmp [lindex $cov [expr 3*$i + 2]]
	    lappend coords $tmp
	    unset tmp
	}
    }


    # Place the beads in preliminary positions on the unit sphere
    for { set i 1 } { $i <= $natomscov } { incr i } {
	set tmp [lindex $coords [expr $i -1]]
	set partnum [lindex $mol  $i ]
	set parttype [lindex $typeinfo 2 [expr $i -1]]
#	puts "$i $partnum $parttype"
	part $partnum pos [lindex $tmp 0] [lindex $tmp 1] [ lindex $tmp 2] type $parttype
    }

    # Make a list of unique particle types in the sphere
    set atomtypelist [::mbtools::utils::uniquelist [lrange $atomtypes 0 [expr $natomscov -1] ] ]

    # Calculate the minimum distance between beads on the sphere which
    # should be roughly equal to the bond length
    set mdist [analyze mindist  $atomtypelist $atomtypelist ]

    # Based on a desired value of mdist equal to 1.0 find the required radius
    set radius [expr 1.0/(1.0*$mdist)]
    mmsg::send [namespace current] "creating hollow sphere with radius $radius and $natomscov beads"


    # Find the total number of each atom type in this hollowsphere
    set sproportions [::mbtools::utils::calc_proportions [lindex $typeinfo 2]]
    set nbeads1 [lindex $sproportions 0 1]
    set nbeads2 [lindex $sproportions 1 1]
    set janus 0
    # If there are exactly 3 types of atoms we assume we have a "Janus Colloid"
    if { [llength $sproportions] == 3 && [expr $nbeads1 + $nbeads2] == $natomscov } {
	# Find the height cutoff (relative to center) for the bead type boundary
	set totalarea [expr 4*3.14159*$radius*$radius]
	set fraction1 [expr $nbeads1/(1.0*$natomscov)]
	set fraction2 [expr $nbeads2/(1.0*$natomscov)]
	set areaf2 [expr $totalarea*$fraction2]
	set heightcut [expr $radius - ($areaf2/(2.0*3.141592*$radius))]
	set janus 1
    } 

    foreach point $coords {
	# Shift the center of the sphere to pos and expand to radius
	for { set x 0 } { $x < 3 } { incr x } {
	    lset point $x [expr [lindex $point $x]*$radius + [lindex $pos $x]]
	}
	lappend scaledcoords $point
    }

   
    set beadcount1 0
    set beadcount2 0
    set tries 0
    set hchange 0.1
    while { $beadcount2 != $nbeads2 && $tries < 100 } {
	set beadcount1 0
	set beadcount2 0
	# Now replace the beads in their new positions
	for { set i 1 } { $i <= $natomscov } { incr i } {
	    set tmp [lindex $scaledcoords [expr $i -1]]
	    set partnum [lindex $mol  $i ]
	    set parttype [lindex $typeinfo 2 [expr $i -1]]
	    # But if we have a janus colloid then use a different technique 
	    # to allocate the bead types
	    if { $janus } {
		# First find the projection of the point onto the unit orient vector
		set proj [::mbtools::utils::dot_product  $tmp  $orient ]
		if { $proj > $heightcut } {
		    set parttype [lindex $sproportions 1 0 ]
		    incr beadcount2
		} else {
		    set parttype [lindex $sproportions 0 0 ]
		    incr beadcount1
		}
	    } else {
		set beadcount2 $nbeads2
	    }

	    part $partnum pos [lindex $tmp 0] [lindex $tmp 1] [ lindex $tmp 2] type $parttype
	    part $partnum fix 1 1 1
	}

	if { $hchange < 0 } {
	    if { $beadcount2 < $nbeads2 } {
		set hchange [expr $hchange*1.5]
	    } else {
		set hchange [expr $hchange*(-0.5)]
	    }
	} else {
	    if { $beadcount2 > $nbeads2 } {
		set hchange [expr $hchange*1.5]
	    } else {
		set hchange [expr $hchange*(-0.5)]
	    }
	}

	if { $janus } {
	    set heightcut [expr $heightcut + $hchange]
	}
#	puts "$beadcount1 $beadcount2 $heightcut"
	incr tries
    }

    if { ($beadcount1 != $nbeads1 || $beadcount2 != $nbeads2) && $janus } {
	mmsg::warn [namespace current] "deviation from intended janus colloid composition $nbeads1 $nbeads2 desired but $beadcount1 $beadcount2 obtained"
    }

    # Now check that the value of mdist is 1.0
    set mdist [analyze mindist  $atomtypelist $atomtypelist ]
    set tol 0.0001
    if { [expr ($mdist - 1.0)*($mdist - 1.0) ] > $tol } {
	mmsg::err [namespace current] "min bond length on sphere is $mdist but it should be 1.0 with a tolerance of $tol"
    }

    
    set mdisttol 0.5

    # Figure out which beads are nearest neigbours
    for { set i 1 } { $i <= $natomscov } { incr i } {
	set partnum [lindex $mol  $i ]
	lappend nbrs [analyze nbhood $partnum [expr $mdisttol + $mdist]]
    }

    # Bond the nearest neighbours
    for { set i 1 } { $i <= $natomscov } { incr i } {
	set partnum [lindex $mol  $i ]
	set nblist [lindex $nbrs [expr $i - 1]]
	foreach atom $nblist {	    
	    if { $atom > $partnum } {
#		Only place the bond if we haven't already done so.
		part $partnum bond $bonds $atom
	    }
	}
    }


    # ------- Now that we have a nice cage we need to fill it with soft balls ----- #
    for { set i $natomscov } { $i < [llength $atomtypes  ] } { incr i } {

	set atomid [lindex $mol [expr $i + 1]]
	set atomtype [lindex $atomtypes $i]

	set maxtries 10000
	set tries 0
	set bpos { 0 0 0 }
	set isallowed 0
	set rbuff 1.0


	while { !$isallowed  } {

	    if {  ($tries > $maxtries) } {
		mmsg::err [namespace current] "could not place molecule exceeded max number of tries"
	    }


	    # First we choose a random point in space within a cube
	    # centered at the center of the sphere
	    lset bpos 0 [expr $radius*2*([t_random] -0.5) + [lindex $pos 0]]
	    lset bpos 1 [expr $radius*2*([t_random] -0.5) + [lindex $pos 1]]
	    lset bpos 2 [expr $radius*2*([t_random] -0.5) + [lindex $pos 2]]
	    # Check to see if the random point is within our sphere and a tolerance
	    if { [mbtools::utils::distance $bpos $pos] < [expr $radius - $rbuff]  } {

		# Now check to see that the bead is within a mindist
		# of other beads in the colloid because these will not
		# be allowed to warm up normally
		part $atomid pos [lindex $bpos 0] [lindex $bpos 1] [lindex $bpos 2] type $atomtype
		if { [analyze mindist $atomtype $atomtype] > 0.6 } {
		    set isallowed 1
		} else {
		    mmsg::debug [namespace current] "tried $tries to place atomid: $atomid"
		}
	    }
	
	    incr tries
	
	}

	part $atomid pos [lindex $bpos 0] [lindex $bpos 1] [lindex $bpos 2] type $atomtype 
	part $atomid fix 1 1 1
    }


    return

}

