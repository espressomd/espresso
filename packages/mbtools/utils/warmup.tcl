   
namespace eval ::mbtools::utils {
    namespace export warmup  
    namespace export nptwarmup
}

#----------------------------------------------------------#
# ::mbtools::utils::warmup--
#
# basic warmup routine
#
# Can be used in a variety of ways depending on the parameters and
# options that are set. 
#
# Options:
#
# mindist: Set this option to a minimum particle distance requirement.
# Routine will use the mindist function of expresso (very very slow )
# to determine if this criterion is satisfied and if it is the warmup
# will be terminated.
#
# imd: Use imd to watch the warmup on the fly 
#
# capincr: This is normally calculated from capgoal and startcap.  If
# it is set then it will override the calculated value
#
proc ::mbtools::utils::warmup { steps times args } {
    variable warmcfg
    ::mmsg::send [namespace current] "warming up $times times $steps timesteps "
    set options {
	{mindist.arg     0   minimum distance criterion to terminate warmup }
	{cfgs.arg -1   write out a record of the warmup at given intervals}
	{outputdir.arg  "./" place to write warmup configs}
	{vmdflag.arg "offline" vmd settings}
	{startcap.arg     5    initial forcecap}	
	{capgoal.arg      1000 forcecap goal}
    }
    set usage "Usage: warmup steps times \[mindist:imd:startcap:capincr:capgoal:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    #Make a sanity check
    if { $steps == 0 || $times == 0 } {
	::mmsg::warn [namespace current] "warmup steps are zero"
	return
    }

    # Work out the cap increment if it is not set
    set capincr [expr ($params(capgoal) - $params(startcap))/($times*1.0)]


    # Set the initial forcecap
    set cap $params(startcap)


    for { set i 0 } { $i < $times } { incr i } {
	# Check the mindist criterion if necessary
	if { $params(mindist) } {
	    set act_min_dist [analyze mindist]
	    if { $act_min_dist < $params(mindist) } {
		break
	    }
	}

	# Write out configuration files and pdb files
	if { $i%$params(cfgs)==0 && ($params(cfgs) > 0 ) } {
	    polyBlockWrite "$params(outputdir)/warm.[format %04d $warmcfg].gz" {time box_l npt_p_diff } {id pos type v f molecule} 
	    mmsg::send [namespace current] "wrote file $params(outputdir)/warm.[format %04d $warmcfg].gz " 

	    flush stdout

	    if { $params(vmdflag) == "offline" } {
		writepdbfoldtopo "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  
	    }
	    incr warmcfg
	}

	# Set the new forcecap into espresso and integrate
	inter tabforcecap $cap
	inter ljforcecap $cap
	integrate $steps
	set cap [expr $cap + $capincr ]
	::mmsg::send [namespace current]  "run $i of $times at time=[setmd time] (cap=$cap) " 

	flush stdout
	
    }
    
    # Note that if we don't reach our desired cap as when capincr is
    # set manually then we simply set it so at the end.
    if { [expr $cap - $capincr] < $params(capgoal) } {
	::mmsg::send [namespace current] "setting final forcecap to $params(capgoal)"
	inter tabforcecap $params(capgoal)
	inter ljforcecap $params(capgoal)
    }
}


#----------------------------------------------------------#
# ::mbtools::utils::nptwarmup--
#
# Routine for decreasing the box mass and friction and allowing box to reach
#
# Can be used in a variety of ways depending on the parameters and
# options that are set. 
#
#
# Note: Assumes a 2d npt ie constant tension
# Options:
#
# 
#
# 
proc ::mbtools::utils::nptwarmup { steps times iparms fparms args } {
    variable warmcfg
    ::mmsg::send [namespace current] "nptwarmup $times times $steps timesteps "
    set options {
	{cfgs.arg -1   write out a record of the warmup at given intervals}
	{outputdir.arg  "./" place to write warmup configs}
	{vmdflag.arg "offline" vmd settings}
    }
    set usage "Usage: warmup steps times \[mindist:imd:startcap:capincr:capgoal:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    #Make a sanity check
    if { $steps == 0 || $times == 0 } {
	::mmsg::warn [namespace current] "npt warmup steps are zero"
	return
    }

    for { set p 0 } { $p < [llength  $iparms] } { incr p } {
	lappend parmincs [expr ([lindex $fparms $p ] - [lindex $iparms $p] )/(1.0*$times-1.0)]

    }

    set systemtemp [setmd temperature ]

    for { set i 0 } { $i < $times } { incr i } {

	# Write out configuration files and pdb files
	if { $i%$params(cfgs)==0 && ($params(cfgs) > 0 ) } {
	    polyBlockWrite "$params(outputdir)/warm.[format %04d $warmcfg].gz" {time box_l npt_p_diff } {id pos type v f molecule} 
	    mmsg::send [namespace current] "wrote file $params(outputdir)/warm.[format %04d $warmcfg].gz " 

	    flush stdout

	    if { $params(vmdflag) == "offline" } {
		writepdbfoldtopo "$params(outputdir)/warm.vmd[format %04d $warmcfg].pdb"  
	    }
	    incr warmcfg
	}

	set p_ext [expr [lindex $iparms 0] + $i*[lindex $parmincs 0] ]
	set piston_mass [expr [lindex $iparms 1] + $i*[lindex $parmincs 1] ]	
	set gamma_0 [expr [lindex $iparms 2] + $i*[lindex $parmincs 2] ]
	set gamma_v [expr [lindex $iparms 3] + $i*[lindex $parmincs 3] ]

	integrate set npt_isotropic $p_ext $piston_mass 1 1 0 
	thermostat set npt_isotropic $systemtemp $gamma_0 $gamma_v 
	mmsg::send [namespace current] "npt integrator has been set $p_ext $piston_mass $systemtemp $gamma_0 $gamma_v "

	integrate $steps

	::mmsg::send [namespace current]  "run $i of $times at time=[setmd time] " 

	flush stdout
	
    }
    
}



## -------- Outdated stuff -------------------- ## 


proc ::mbtools::utils::free_warmup { args } {
    ::mmsg::send [namespace current] "warming up .. " nonewline
    set options {
	{bondl.arg     1.0   bond length between atoms  }
	{uniform.arg    1    use uniform lipid placement }
	{nhb.arg        1    number of head beads per lipid}	
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:uniform]"
    array set params [::cmdline::getoptions args $options $usage]
    
    


    set steps 1000
    set warm_n_times 10
    set imd_output "off"
    
    set count 0
    set skip 0
    foreach a $args {
	set next [expr $count + 1]
	if { $skip == 0} {
	    if { $a == "steps" } { 
		set steps  [lindex $args $next] 
	    } elseif { $a == "times" } { 
		set times [lindex $args $next] 
	    } elseif { $a == "criterion" } { 
		set criterion  [lindex $args $next] 
	    } elseif { $a == "imd_output" } { 
		set imd_output  [lindex $args $next] 
	    }  else {
		return -code error -errorinfo "$a is not an argument of warmup";
	    }
	    set skip 1
	} else { set skip 0 }
	set count [expr $count + 1]
    }
    
    set i 0
    ::mmsg::send [namespace current]  "warming up: " nonewline
    while { $i < $times } {
	::mmsg::send [namespace current]  ". " nonewline
	flush stdout
	integrate $steps
	
	# One day make a proper test
	#	if { [lindex $criterion 0] == "pressure" } {
	#	    lappend idealp [lindex [analyze pressure ideal]  1] 
	#	    
	#	}
	
	
	# Visualization
	if { $imd_output=="interactive" } { imd positions -fold_chains }
	incr i
    }
    
}