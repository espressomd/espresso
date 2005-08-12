# system_generation --
#
# Routines in this package provide a variety of premade methods for
# setting up typical systems such a vesicles toroids and bilayers
# Author: Ira
# 

package require mmsg 0.1.0
package require mathutils 0.1.0
package provide system_generation 1.0.0

# Create the empty namespace to which we shall add our routines
namespace eval ::system_generation {
    # Global routine for setting up the system
    namespace export setup_system

    # Global variables for setup routines
    variable moltypeskey  
    variable userfixedparts

    variable middlebead

    variable icovermagicnums { 72 92 122 132 162 192 212 252 272 282 312 362 372 392 432 482 492 492 522 572 612 632 642 672 732 752 762 792 812 842 912 912 932 972 1002 1032 1082 1092 1112 1122 1172 1212 1242 1272 1292 1332 1332 1392 1442 1472 1472 1482 1512 1562 1572 1632 1692 1692 1712 1722 1752 1812 1832 1892 1922 1932 1962 1962 1992 2012 2082 2112 2172 2172 2192 2232 2252 2282 2292 2372 2412 2432 2442 2472 2472 2522 2562 2592 2592 2682 2712 2732 2732 2772 2792 2832 2892 2912 2922 3002 3012 3012 3042 3072 3092 3132 3162 3242 3252 3272 3312 3332 3362 3372 3432 3432 3492 3512 3612 3612 3632 3642 3642 3672 3722 3732 3792 3812 3872 3882 3972 3992 3992 4002 4032 4032 4092 4122 4172 4212 4272 4272 4322 4332 4362 4392 4412 4412 4442 4482 4532 4572 4632 4682 4692 4692 4712 4752 4812 4812 4842 4872 4892 4962 4992 5072 5072 5082 5112 5112 5132 5162 5232 5252 5292 5322 5322 5412 5432 5472 5492 5532 5532 5562 5592 5592 5672 5712 5762 5772 5792 5882 5882 5892 5892 5922 5972 6012 6032 6042 6072 6132 6192 6242 6252 6282 6312 6332 6372 6372 6372 6432 6512 6512 6522 6572 6612 6692 6732 6752 6762 6762 6792 6792 6842 6872 6882 6912 7002 7682 8192 8192 8672 9722 10832 12002 13232 14522 15872 17282 18752 20282 21872 23522 25232 27002 28832 30722 32672 34682 36752 38882 41072 43322 45632 48002 50432 52922 55472 58082 60752 63482 66272 69122 72032 75002 78032 }


    # Read in all the routines
    source [file join [file dirname [info script]] bilayer.tcl]
    source [file join [file dirname [info script]] random.tcl]
    source [file join [file dirname [info script]] create_torus.tcl]
    source [file join [file dirname [info script]] sphere.tcl]
    source [file join [file dirname [info script]] cylinder.tcl]
    source [file join [file dirname [info script]] constraint.tcl]
    source [file join [file dirname [info script]] topologies.tcl]
    source [file join [file dirname [info script]] utils.tcl]
    source [file join [file dirname [info script]] place.tcl]
    source [file join [file dirname [info script]] readfromfile.tcl]
    



}

# ::system_generation::setup_system -- 
#
# A large scale wrapper routine for all of the system setup commands
# in this module.
# 
# This routine should allow flexible use of multiple geometrical
# structures and topologies in a simulation.  A complete topology and
# geometry specification is provided for each object (eg sphere, flat
# bilayer etc) in the system and these are grouped into the list
# <system_specs>.  setup_system loops through each of these objects
# setting each up independently and then combining topologies at the
# end.
#
# Arguments:
# 
# system_specs: A list of specified system objects.  Each object is
# itself a list of the form < geometry n_lipidslist beads_per_mol >
# where geometry may be any allowable geometry (eg flat sphere
# etc). <n_lipidslist> is a list containing the number of molecules of
# each of the molecule types in this object and beads_per_mol
# specifies the number of particles in each of the molecule types.
#
# 
# 
#
#
proc ::system_generation::setup_system { system_specs setbox_l moltypes } {
    ::mmsg::send [namespace current] "setting up system "

    # The molecule types spec should be globally accessible
    variable moltypeskey
    variable userfixedparts
    set moltypeskey $moltypes

    # Starting value for particle ids
    set startp 0
    # Current particle id
    set currpid $startp


    foreach spec $system_specs {
	# Flags for tracking input settings
	set geometryset 0
	set n_molslistset 0
	set n_molslist 0
	set notopo 0
	foreach item $spec {
	    switch [lindex $item 0] {
		"geometry" {
		    set geometry [lindex $item 1]
		    set geometryset 1
		}
		"n_molslist" {
		    set n_molslist [lindex $item 1]
		    set n_molslistset 1
		}
		"default" {
		    ::mmsg::warn [namespace current] "unknown item [lindex $item 0] in system spec. allowed values are: \n geometry \n n_molslist  "
		}
	    }
	}



	if { !$geometryset } {
	    mmsg::err [namespace current] "geometry not specified"
	}
	if { !$n_molslistset } {
	    mmsg::err [namespace current] "n_molslist not specified for [lindex $geometry 0]"
	}

	# Generate a topology from an a list of the number and size of
	# each molecule
	foreach mtp $n_molslist {
	    set thismoltypeid [lindex $mtp 0]
	    set nmols [lindex $mtp 1]
	    set tpspec [matchtype [lindex $mtp 0]]
	    set nbeads_mol [llength [lindex $tpspec 2] ]
	    # Create the topology for this lipid type
	    set topo [::system_generation::create_simple_topo $nmols $nbeads_mol -moltype  $thismoltypeid -startpart $currpid ]		
		
	    # Just in case zero molecules were specified we need
	    # to check if topo was actually created at all by the
	    # last command
	    if { $topo == 0 } {
		::mmsg::err [namespace current] "no topo created for molecule type $thismoltypeid"
	    } else {
		lappend topolist $topo
		set currpid [expr [maxpartid $topo ] + 1]
	    }

	}

	# Join all of the previously made topologies
	set first 1
	foreach topo $topolist {
	    if { $first } { 
		set topology $topo 
		set first 0
	    } else {
		set topology [::system_generation::join_topos $topology $topo]
	    }
	    
	}
	unset topolist


	# Now wrap the topology onto a specified geometry and perform any
	# other geometry specific tasks
	switch  [lindex $geometry 0]  {
	    "flat" {

		::system_generation::create_bilayer $topology $setbox_l 

		# Fix the z positions for warmup
		::mmsg::send [namespace current] "fixing z positions of bilayer for warmup" 
		for {set i [minpartid $topology] } { $i <  [setmd n_part] } {incr i} {
		    set fixvalue [part $i print fix]
		    if { [lindex $fixvalue 0] == 0 && [lindex $fixvalue 1] == 0 && [lindex $fixvalue 2] == 0 }  {
			part [expr $i] fix 0 0 1
		    } else {
			lappend userfixedparts $i
		    }
		}  

		# Check particle consistency
		if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		    mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}
	    }
	    "sphere" {
		if { [llength $geometry] > 1 } {
		    set center [lindex $geometry 1]
		} else {
		    set center [list 0.0 0.0 0.0]
		}
		set topology [::system_generation::shuffle_topo $topology ]
		::system_generation::create_sphere $topology $setbox_l -center $center 

		# Check particle consistency
		if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		    mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}
		
	    }
	    "cylinder" {
		if { [llength $geometry] > 1 } {
		    set center [lindex $geometry 1]
		} else {
		    set center [list 0.0 0.0 0.0]
		}
		set topology [::system_generation::shuffle_topo $topology ]
		::system_generation::create_cylinder $topology $setbox_l -center $center

		# Check particle consistency
		if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		    mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}

	    }
	     "random" {

		 if { [llength $geometry] > 1 } {
		     set extraargs [lrange $geometry 1 end]
		 } else { set extraargs 0} 



		 ::system_generation::create_random_fluid $topology $setbox_l  -exclude [concat $extraargs] 

		# Check particle consistency
		 if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		     mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}
	    }
	    "singlemol" {
		if { [llength $geometry ] > 2 } {
		    set orient [lindex $geometry 2]
		} else {
		    set orient { 0 0 1 }
		}
		set center [lindex $geometry 1]
		set mol [lindex $topology 0]
		placemol $mol $center -orient $orient

		# Retrieve the molecule information for this molecule type
		set typekey [matchtype [lindex $mol 0]]

		# Check particle consistency
		if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		    mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}

		# Special condition for constraints is that they
		# should not contribute to the topology
		if  { [lindex $typekey 1] == "sphericalconstraint" } {
		    set notopo 1
		    unset topology
		}
	    }
	    "readfile" {
		set filename [lindex $geometry 1]
		set topofile [lindex $geometry 2]
		set topology [::system_generation::readfromfile $topology $setbox_l $filename $topofile -ignore { f } ]

		# Check particle consistency
		if { [setmd n_part] != [expr [maxpartid $topology] + 1] } {
		    mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [maxpartid $topology] +1] were specified in topology "
		}
	    }
	    "default" {
		mmsg::err [namespace current] "no setup method known for [lindex $geometry 1]"
	    }
	    
	}

	if {!$notopo} {
	    lappend topologieslist $topology
	}
	
    }

    # Join all of the previously made topologies
    set first 1
    foreach topo $topologieslist {
	if { $first } { 
	    set topology $topo 
	    set first 0
	} else {
	    set topology [::system_generation::join_topos $topology $topo]
	}
	
    }

    set topology [::system_generation::sort_topo $topology]

    return $topology

}



proc ::system_generation::get_userfixedparts {  } {
    variable userfixedparts

    if { [catch { set dum $userfixedparts } ] } {
	::mmsg::warn [namespace current] "no user fixed particles defined"
    }
    return $userfixedparts
}

proc ::system_generation::get_middlebead {  } {
    variable middlebead
    return $middlebead
}
