#
#
# Routines for placing molecules
#



namespace eval system_generation {


}



# ::system_generation::placemol-- 
#
# general routine for placing molecules
#
#
proc ::system_generation::placemol { mol pos args } { 
    set options {
	{bondl.arg     1.0   "bond length between atoms"  }
	{orient.arg  { 0 0 1 } "orientation vector for the mol " }
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:orient]"
    array set params [::cmdline::getoptions args $options $usage]
    

    # Retrieve the molecule information for this molecule type	
    set typekey [matchtype [lindex $mol 0]]

    # Place the molecule depending on type
    switch [lindex $typekey 1] {
	"lipid" {
	    place_lipid_linear $mol $params(orient) $pos -bondl $params(bondl) 
	}
	"spanlipid" {
	    place_lipid_linear $mol $params(orient) $pos -bondl $params(bondl) -midpos 
	}
	"hollowsphere" {
	    place_hollowsphere $mol $pos
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

# ::system_generation::matchtype-- 
#
# Search the molecule type list for a type number and return the type specification
#
#
proc ::system_generation::matchtype { mol } {

    variable moltypeskey

    foreach key $moltypeskey {
	if { [lindex $key 0] == [lindex $mol 0] } {
	    return $key
	}
    }
    mmsg::err [namespace current] "could not find a matching key to moltype [lindex $mol 0]"
}



# ::system_generation::place_lipid_linear-- 
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
proc ::system_generation::place_lipid_linear { mol orient pos args } {
    set options {
	{bondl.arg     1.0   bond length between atoms  }
	{midpos  "the value of pos corresponds to the middle bead not the tail"}
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:uniform]"
    array set params [::cmdline::getoptions args $options $usage]



    # Ensure that the orientation vector is normalized
    set orient [::mathutils::normalize $orient]


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
	# placement from the tail we need to use the following
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


# ::system_generation::place_sphericalconstraint-- 
#
# Place a spherical constraint
#
# Arguments:
#
# mol: Information on the particle type for the constraint
# pos: The center of the sphere
#
proc ::system_generation::place_sphericalconstraint { mol pos args } {
    set options {
    }
    set usage "Usage: place_lipid_linear topo boxl \[mol:pos:]"
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
proc ::system_generation::place_protein { mol pos orient args } {

    variable middlebead

    set moltype [lindex $mol 0]
    set typeinfo [matchtype $moltype]

    set atomtypes [lindex $typeinfo 2]
    set bonds [lindex $typeinfo 3]

  # Ensure that the orientation vector is normalized
    set orient [::mathutils::normalize $orient]

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




# ::system_generation::place_hollowsphere-- 
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
proc ::system_generation::place_hollowsphere { mol pos args } {
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

#    puts $typeinfo




    # Place the beads in preliminary positions on the unit sphere
    for { set i 1 } { $i <= $natomscov } { incr i } {
	set tmp [lindex $coords [expr $i -1]]
	set partnum [lindex $mol  $i ]
	set parttype [lindex $typeinfo 2 [expr $i -1]]
#	puts "$i $partnum $parttype"
	part $partnum pos [lindex $tmp 0] [lindex $tmp 1] [ lindex $tmp 2] type $parttype
    }

    # Make a list of unique particle types in the sphere
    set atomtypelist [uniquelist [lrange $atomtypes 0 [expr $natomscov -1] ] ]

    # Calculate the minimum distance between beads on the sphere which
    # should be roughly equal to the bond length
    set mdist [analyze mindist  $atomtypelist $atomtypelist ]

    # Based on a desired value of mdist equal to 1.0 find the required radius
    set radius [expr 1.0/(1.0*$mdist)]
    mmsg::send [namespace current] "creating hollow sphere with radius $radius and $natomscov beads"


    # Find the total number of each atom type in this hollowsphere
    set sproportions [::mathutils::calc_proportions [lindex $typeinfo 2]]
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
		if { [expr [lindex $tmp 2] - [lindex $pos 2] ] > $heightcut } {
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

	set heightcut [expr $heightcut + $hchange]
#	puts "$beadcount1 $beadcount2 $heightcut"
	incr tries
    }

    if { $beadcount1 != $nbeads1 || $beadcount2 != $nbeads2 } {
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
	    if { [mathutils::distance $bpos $pos] < [expr $radius - $rbuff]  } {

		# Now check to see that the bead is within a mindist
		# of other beads in the colloid because these will not
		# be allowed to warm up normally
		part $atomid pos [lindex $bpos 0] [lindex $bpos 1] [lindex $bpos 2] type $atomtype
		if { [analyze mindist $atomtype $atomtype] > 0.6 } {
		    set isallowed 1
		} else {
		    mmsg::send [namespace current] "tried $tries to place atomid: $atomid"
		}
	    }
	
	    incr tries
	
	}

	part $atomid pos [lindex $bpos 0] [lindex $bpos 1] [lindex $bpos 2] type $atomtype 
	part $atomid fix 1 1 1
    }


    return

}


