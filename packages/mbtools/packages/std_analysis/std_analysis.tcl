
package require mmsg 0.1.0
package provide std_analysis 0.1.0

namespace eval ::std_analysis {

    # Global arguments
    variable area_lipid
    variable iotype 
    variable all_particles
    variable rawmodes
    variable hrange
    variable lhrange
    variable lorange 
    variable beadtypes
    variable nbins
    variable profilenogrid
    variable l_orients_start
    variable suffix
    variable topology

    # Special or hardwired variables
    variable rcatch 1.9

    variable localheightsnbins 100
    variable localorientsnbins 100




    variable switches
    variable this [namespace current]
    variable known_flags " possible flags are: \n cluster_calc \n pik1_calc \n pressure_calc \n box_len_calc \n fluctuation_calc \n energy_calc \n stray_lipids_calc \n orient_order_calc \n flipflop_calc \n density_profile_calc \n localheights_calc \n localorients \n distance_calc \n "

    #File Streams
    variable f_tvspik1
    variable f_tvsp
    variable f_tvsbl
    variable f_tvsflip
    variable f_tvsoop
    variable f_tvsstray
    variable f_tvsen
    variable f_tvsclust
    variable f_densityprof
    variable f_localheights
    variable f_localorients
    variable f_tvsdist
    # Variables to be used for averaging
    variable av_localorients 0
    variable av_localorients_i 0

    variable av_dist  0
    variable av_dist_i 0


    variable av_densities
    variable av_densities_i 0

    variable av_localheights
    variable av_localheights_i 0

    variable av_pow
    variable av_pow_i 0

    variable av_sizehisto 0
    variable av_sizehisto_i 0

    variable av_clust { 0 0 0 0 0 0 0 0 0 0 0 }
    variable av_clust_i 0

    variable av_pik1 { 0 0 0 0 0 0 0 0 0 }
    variable av_pik1_i 0

    variable av_pressure { 0 0 0 0 0 0 }
    variable av_pressure_i 0

    variable av_boxl { 0 0 0 }
    variable av_boxl_i 0

    variable av_flip 0.0
    variable av_flip_i 0

    variable av_oop 0.0
    variable av_oop_i 0
    
    variable av_stray 0
    variable av_stray_i 0

    variable av_components_en 0
    variable av_total_en 0
    variable av_kin_en 0
    variable av_fene_en 0
    variable av_harm_en 0
    variable av_nb_en 0
    variable av_en_i 0


    namespace export do_analysis
    namespace export setup_analysis
    namespace export print_averages
    namespace export flush_streams

}

source [file join [file dirname [info script]] flipflop.tcl]
source [file join [file dirname [info script]] boxl.tcl]
source [file join [file dirname [info script]] clusters.tcl]
source [file join [file dirname [info script]] energy.tcl]
source [file join [file dirname [info script]] pressure.tcl]
source [file join [file dirname [info script]] pik1.tcl]
source [file join [file dirname [info script]] oop.tcl]
source [file join [file dirname [info script]] fluctuations.tcl]
source [file join [file dirname [info script]] stray.tcl]
source [file join [file dirname [info script]] density_profile.tcl]
source [file join [file dirname [info script]] localheights.tcl]
source [file join [file dirname [info script]] localorients.tcl]
source [file join [file dirname [info script]] distance.tcl]

# ::std_analysis::flush_streams --
#
# Flush all of the output stream associated with std_analysis
#
proc ::std_analysis::flush_streams { } {
    variable known_flags
    variable switches
    variable this
    variable f_tvspik1
    variable f_tvsp
    variable f_tvsbl
    variable f_tvsflip
    variable f_tvsoop
    variable f_tvsen
    variable f_tvsclust
    variable f_tvsstray
    variable f_densityprof
    variable f_localheights
    variable f_localorients
    variable f_tvsdist

    for { set i 0 } { $i < [llength $switches ] } { incr i } {
	switch [lindex $switches $i 0] {
	    "fluctuation_calc" {
		# This doesn't get flushed because we overwrite each time
	    }
	    "density_profile_calc" {
		# This doesn't get flushed because we overwrite each time
	    }
	    "localheights_calc" {
		# This doesn't get flushed because we overwrite each time
	    }
	    "localorients_calc" {
		# This doesn't get flushed because we overwrite each time
	    }
	    "cluster_calc" {
		flush $f_tvsclust
	    }
	    "pik1_calc" {
		flush $f_tvspik1
	    }
	    "pressure_calc" {
		flush $f_tvsp
	    }
	    "box_len_calc" {
		flush $f_tvsbl
	    }
	    "flipflop_calc" {
		flush $f_tvsflip
	    }
	    "orient_order_calc" {
		flush $f_tvsoop
	    }
	    "energy_calc" {
		flush $f_tvsen
	    }
	    "stray_lipids_calc" {
		flush $f_tvsstray
	    }
	    "distance_calc" {
               flush $f_tvsdist
	    }
	    "default" {
		mmsg::warn $this "unknown analysis flag [lindex $switches $i 0] $known_flags" 
	    }
	}	    
    }
}

# ::std_analysis::print_averages --
#
# Calculate averages for all analyzed quantities and put them to the
# appropriate file streams
#
proc ::std_analysis::print_averages { } {
    variable this
    mmsg::debug $this "printing averages"
    variable switches
    variable known_flags
    variable topology
    variable localheightsnbins
    variable lhrange
    variable localorientsnbins
    variable lorange

    variable av_densities
    variable av_densities_i 
    variable f_densityprof


    variable av_localheights
    variable av_localheights_i
    variable f_localheights

    variable av_localorients
    variable av_localorients_i
    variable f_localorients


    variable av_pik1_i
    variable av_pik1
    variable f_tvspik1

    variable av_clust
    variable av_clust_i
    variable f_tvsclust

    variable av_pressure
    variable av_pressure_i
    variable f_tvsp

    variable av_boxl_i
    variable av_boxl
    variable f_tvsbl

    variable av_pow_i
    variable av_pow
    
    variable av_flip
    variable av_flip_i
    variable f_tvsflip
    
    variable av_oop
    variable av_oop_i
    variable f_tvsoop
    
    variable av_stray
    variable av_stray_i
    variable f_tvsstray
    
    variable av_dist_i
    variable av_dist
    variable f_tvsdist

    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_fene_en
    variable av_harm_en
    variable av_nb_en
    variable av_en_i
    
    
    variable f_tvsen
    
    variable av_sizehisto
    variable av_sizehisto_i
    variable outputdir
    variable iotype
    variable suffix
    variable beadtypes

    variable hrange
    variable nbins
    
    set binwidth [expr $hrange*2.0/(1.0*$nbins)]
    
    set time [setmd time]
    
    mmsg::debug $this "calculating [llength $switches ] different quantities"
    for { set i 0 } { $i < [llength $switches ] } { incr i } {
	#	    puts [lindex $switches $i 0]
	#	    flush stdout
	switch [lindex $switches $i 0] {
	    "fluctuation_calc" {
		if { [lindex $switches $i 2] && $av_pow_i > 0 } {
		    # Do a power analysis on the intermediate results
		    power_analysis
		}
	    }
	    "density_profile_calc" {
		if { [lindex $switches $i 2] && $av_densities_i > 0 } {
		    
		    set f_densityprof [open "$outputdir/av_zprof$suffix" "w" ]
		    # Write a header to the file
		    puts -nonewline $f_densityprof "\# zheight "
		    foreach bt $beadtypes {
			puts -nonewline $f_densityprof "|| $bt up "
		    }
		    foreach bt $beadtypes {
			puts -nonewline $f_densityprof "|| $bt down "
		    }
		    puts $f_densityprof ""

		    # Write out the density profiles to a file
		    for { set bin 0 } { $bin < [llength $av_densities] } { incr bin } {
			set currbin [lindex $av_densities $bin]
			puts -nonewline $f_densityprof "[expr $bin*$binwidth-($binwidth/2.0)] "
			for { set bt 0 } { $bt < [llength $currbin] } { incr bt } {
			    puts -nonewline $f_densityprof "[expr [lindex $currbin $bt]/(1.0*$av_densities_i)] "
			}
			puts $f_densityprof ""
		    }		    
		    close $f_densityprof
		}
	    }
	    "localheights_calc" {
		if { [lindex $switches $i 2] && $av_localheights_i > 0 } {
		    
		    set f_localheights [open "$outputdir/av_localh$suffix" "w" ]
		    # Write a header to the file
		    puts $f_localheights "\# local height distribution of lipids "
		    puts $f_localheights ""

		   
		    set lhbinwidth [expr $lhrange/($localheightsnbins*1.0)]

		    # Write out the local height distribution to a file
		    for { set bin 0 } { $bin < [llength $av_localheights] } { incr bin } {
			set currbin [lindex $av_localheights $bin]
			puts $f_localheights "[expr $bin*$lhbinwidth - $lhrange/2.0] [expr $currbin/(1.0*[llength $topology])]"
		    }		    
		    close $f_localheights
		}
	    }
	    "localorients_calc" {
		if { [lindex $switches $i 2] && $av_localorients_i > 0 } {
		    
		    set f_localorients [open "$outputdir/av_localo$suffix" "w" ]
		    # Write a header to the file
		    puts $f_localorients "\# local orientation distribution of lipids "
		    puts $f_localorients ""

		   
		    set lobinwidth [expr $lorange/($localorientsnbins*1.0)]

		    # Write out the local orientation distribution to a file
		    for { set bin 0 } { $bin < [llength $av_localorients] } { incr bin } {
			set currbin [lindex $av_localorients $bin]
			puts $f_localorients "[expr $bin*$lobinwidth] [expr $currbin/(1.0*[llength $topology])]"
		    }		    
		    close $f_localorients
		}
	    }
	    "cluster_calc" {
		if { [lindex $switches $i 2] && $av_clust_i > 0 } {
		    puts -nonewline $f_tvsclust "$time "
		    for { set v 0 } { $v < [llength $av_clust] } {incr v} {
			puts -nonewline $f_tvsclust "[expr [lindex $av_clust $v]/($av_clust_i*1.0)] "
		    }
		    puts $f_tvsclust " "
		    
		    set tident [expr int(floor($time))]
		    set tmpfile [open "$outputdir/sizehisto.[format %05d $tident]" w ]
		    for {set v 0 } { $v < [llength $av_sizehisto] } { incr v } {
			puts $tmpfile "[expr $v + 1] [expr ([lindex $av_sizehisto $v]/(1.0*$av_sizehisto_i))]"
		    }
		    close $tmpfile
		    
		} else {
		    mmsg::warn $this "can't print average clusters"
		    flush stdout
		}
	    }
	    "pik1_calc" {
		if { [lindex $switches $i 2] && $av_pik1_i > 0 } {
		    puts -nonewline $f_tvspik1 "$time "
		    for { set v 0 } { $v < [llength $av_pik1] } {incr v} {
			puts -nonewline $f_tvspik1 "[expr [lindex $av_pik1 $v]/($av_pik1_i*1.0)] "
			
		    }
		    puts $f_tvspik1 " "
		} else {
		    mmsg::warn $this "can't print average pik1"
		    flush stdout		       
		}
	    }
	    "pressure_calc" {
		if { [lindex $switches $i 2] && $av_pressure_i > 0 } {
		    puts -nonewline $f_tvsp "$time "
		    for { set v 0 } { $v < [llength $av_pressure] } {incr v} {
			puts -nonewline $f_tvsp "[expr [lindex $av_pressure $v]/($av_pressure_i*1.0)] "
			
		    }
		    puts $f_tvsp " "
		} else {
		    mmsg::warn $this "can't print average pressure"
		    flush stdout		       
		}
	    }
	    "box_len_calc" {
		if { [lindex $switches $i 2] && $av_boxl_i > 0 } {
		    set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
		    set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
		    set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
		    puts $f_tvsbl "$time $avblx $avbly $avblz"
		    #			set av_boxl_i 0
		} else {
		    mmsg::warn $this "can't print average box length"
		    flush stdout
		}
		
	    }

	    "flipflop_calc" {
		puts $f_tvsflip "$time [expr $av_flip/(1.0*$av_flip_i)]"
	    }
	    "orient_order_calc" {
		puts $f_tvsoop "$time [expr $av_oop/(1.0*$av_oop_i)]"
	    }
	    "stray_lipids_calc" {
		puts $f_tvsstray "$time [expr $av_stray/(1.0*$av_stray_i)]"
	    }
	    "energy_calc" {
		puts -nonewline $f_tvsen "$time [expr $av_total_en/(1.0*$av_en_i)] [expr $av_fene_en/(1.0*$av_en_i)] [expr $av_kin_en/(1.0*$av_en_i)] [expr $av_harm_en/(1.0*$av_en_i)] [expr $av_nb_en/(1.0*$av_en_i)]"
		foreach comp $av_components_en {
		    puts -nonewline $f_tvsen " [expr $comp/(1.0*$av_en_i)]"
		}
		puts $f_tvsen ""
	    }
	    "distance_calc" {
               if { [lindex $switches $i 2] && $av_dist_i > 0 } {
                   set avdistx [expr [lindex $av_dist 0]/($av_dist_i*1.0)]

                   puts $f_tvsdist "$time $avdistx"
                   #                   set av_boxl_i 0
               } else {
                   mmsg::warn $this "can't print average distance"
                   flush stdout
               }
	    }
	    "default" {
		mmsg::warn $this "unknown analysis flag [lindex $switches $i 0] $known_flags" 
	    }
	}
    }
    reset_averages
    #	puts "done"
    flush stdout
}


# ::std_analysis::reset_averages --
#
# Reset all of the average storage variables and counters to zero
#   
#
# Note: Power analysis and densityprofiles are not reset since they
# generally require averages over the entire simulation. Flip-flop is
# also not reset since it should exponentially decay with time and is
# calculated from the entire simulation run.
#
proc ::std_analysis::reset_averages { } {
    variable av_sizehisto
    variable av_sizehisto_i
    variable known_flags

    variable av_clust
    variable av_clust_i
    
    variable av_pik1_i
    variable av_pik1
    
    variable av_pressure
    variable av_pressure_i
    
    variable av_boxl_i
    variable av_boxl
    
    variable av_flip
    variable av_flip_i
    
    variable av_oop
    variable av_oop_i
    
    variable av_stray
    variable av_stray_i
    
    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_fene_en
    variable av_harm_en
    variable av_nb_en
    variable av_en_i
    
    variable av_dist_i
    variable av_dist

    set av_dist 0.0
    set av_dist_i 0

    set av_sizehisto 0
    set av_sizehisto_i 0
    
    
    set av_clust {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_clust_i 0
    
    set av_boxl { 0.0 0.0 0.0 }
    set av_boxl_i 0
    
    set av_pik1 { 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_pik1_i 0
    
    set av_pressure {0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_pressure_i 0
    
    set av_flip 0
    set av_flip_i 0
    
    set av_oop 0
    set av_oop_i 0
    
    set av_stray 0
    set av_stray_i 0
    
    for {set i 0 } {$i < [llength $av_components_en] } { incr i } {
	lset av_components_en $i 0.0
    }
    set av_total_en 0
    set av_kin_en 0
    set av_fene_en 0
    set av_harm_en 0
    set av_nb_en 0
    set av_en_i 0
    

    
}

# ::std_analysis::do_analysis --
#
# This is the routine that is typically called during the simulation
# after each integrate command.  It simply calls all of the
# appropriate analysis routines and stores the values in average
# storage variables
#
proc ::std_analysis::do_analysis { } {
    variable switches
    variable this
    variable known_flags

    for { set i 0 } { $i < [llength $switches ] } { incr i } {
	switch [lindex $switches $i 0] {
	    "density_profile_calc" {
		analyze_density_profile [lindex $switches $i 1]
	    }
	    "localheights_calc" {
		analyze_localheights [lindex $switches $i 1]
	    }
	    "localorients_calc" {
		analyze_localorients  [lindex $switches $i 1]
	    }
	    "cluster_calc" {
		analyze_clusters [lindex $switches $i 1]
	    }
	    "pik1_calc" {
		analyze_pik1 [lindex $switches $i 1]
	    }
	    "pressure_calc" {
		analyze_pressure [lindex $switches $i 1]
	    }
	    "quick_pressure_calc" {
		analyze_quick_pressure [lindex $switches $i 1]
	    }
	    "box_len_calc" {
		analyze_box_len [lindex $switches $i 1]
	    }
	    "fluctuation_calc" {
		analyze_fluctuations [lindex $switches $i 1]
	    }
	    "flipflop_calc" {
		analyze_flipflop [lindex $switches $i 1]
	    }
	    "orient_order_calc" {
		analyze_oop [lindex $switches $i 1]
	    }
	    "stray_lipids_calc" {
		analyze_stray [lindex $switches $i 1]
	    }
	    "energy_calc" {
		analyze_energy [lindex $switches $i 1]
	    }
	    "distance_calc" {
               analyze_distance [lindex $switches $i 1]
	    }
	    "default" {
		mmsg::warn $this "unknown analysis flag [lindex $switches $i 0] $known_flags" 
	    }
	}	    
    }
    mmsg::debug $this "done"
    flush stdout
}

# ::std_analysis::setup_analysis --
#
# This routine should be called once and only once at the beginning of
# the simulation in order to setup all the appropriate variables that
# will later be used when do_analysis is called.
#
# Arguments:
#
#        switchesin: This argument should consist of a list of
#                    switches that determine which quantites should be
#                    analyzed and whether results should be printed to
#                    stdout during calculation. To see what this list
#                    looks like just look at the example
#                    parameters.tcl file 
#
#
#       mgrid:       Size of the grid used for calculation of a height grid.
#                    This is used for flip-flop fluctation and stray
#                    lipid routines.
#
#       outputdir:   Output directory
#
#       straycutoff: The Cutoff distance beyond which a lipid is
#                    considered to be a stray
#  
#
#
#        hrange:    The range relative to bilayer midplane over which to
#                   calculate the density profiles
#
#        nbins:     The number of bins used for density profile analysis
#
#       beadtypes:  A list of bead types to be used for density profile analysis
#
#       a_lipid:    The area per lipid
#
#       suffix: The suffix to use for output files.  
#
#       iotype: This parameter set the iotype of the files that will
#                    be opened for writing output.  Set this to "a" if
#                    you want to append to old data or to "w" if you
#                    want to overwrite.
#
proc ::std_analysis::setup_analysis { switchesin topo args } {

    variable suffix
    variable known_flags
    variable iotype
    variable switches
    variable n_particles
    variable mgrid
    variable outputdir
    variable stray_cut_off
    variable all_particles
    variable lhrange
    variable lorange
    variable f_tvspik1
    variable f_tvsp
    variable f_tvsbl
    variable f_tvsflip
    variable f_tvsoop
    variable f_tvsstray
    variable f_tvsen
    variable f_tvsclust
    variable f_tvsdist
    variable av_pow
    variable av_localheights
    variable av_localorients
    variable av_densities
    variable av_densities_i
    variable l_orients_start
    variable topology
    variable area_lipid
    variable hrange
    variable nbins
    variable beadtypes
    variable this
    variable profilenogrid
    variable localheightsnbins
    variable localorientsnbins
    variable av_components_en

    mmsg::send $this "setting up analysis"

    set options {
	{mgrid.arg      8    "set the size of the grid for heightfunction calculations" }
	{straycutoff.arg      4    "stray distance from bilayer " }
	{outputdir.arg      "./"    "name of output directory " }
	{alipid.arg  1.29 "area per lipid" }
	{suffix.arg "tmp" "suffix to be used for outputfiles" }
	{iotype.arg "a" "the method with which to open existing analysis files"}
	{nbins.arg "50" "Number of bins used for density profile analysis" }
	{hrange.arg "5" "Range over which to calculate density profile"}
	{beadtypes.arg "0" "Identity of beads to use for density profile analysis" }
	{profilenogrid.arg "0" "Whether to not use the height grid for density profile analysis" }
	{lhrange.arg "4" "Range for localheights calculation" }
	{lorange.arg "1.0" "Range for localorients calculation" }
    }
    set usage "Usage: setup_analysis gridm:straycutoff:outputdir:alipid:suffix:iotype:nbins:hrange:beadtypes:profilenogrid "
    array set params [::cmdline::getoptions args $options $usage]


    # This checks for a common problem in the formatting of the $switchesin
    # and attempts for construct a list called switches that is
    # correct

    if { [llength [lindex [lindex $switchesin 0] 0 ] ] == 3 } {
	set switches ""
	for { set i 0 } { $i < [ llength [lindex $switchesin 0] ] } {incr i} {
	    lappend switches [lindex [lindex $switchesin 0] $i ]
	}
    } else {
	set switches $switchesin
    }

    # Calculate the total number of particles from topology
    set n_particles 0
    foreach mol $topo {
	set n_particles [expr $n_particles + [llength $mol] -1]
    }

    set mgrid $params(mgrid)
    set outputdir $params(outputdir)
    set stray_cut_off $params(straycutoff)
    set topology $topo
    set iotype $params(iotype)
    set suffix "_$params(suffix)"
    set area_lipid $params(alipid)
    set hrange $params(hrange)
    set nbins $params(nbins)
    set beadtypes $params(beadtypes)
    set profilenogrid $params(profilenogrid)
    set lhrange $params(lhrange)
    set lorange $params(lorange)

    for { set i 0 } { $i < [llength $switches ] } { incr i } {
	mmsg::debug $this "switch = [lindex $switches $i 0]"
	switch [lindex $switches $i 0] {
	    "cluster_calc" {
		mmsg::debug $this "opening $outputdir/time_vs_clust$suffix "		    

		if { [file exists "$outputdir/time_vs_clust$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}

		set f_tvsclust [open "$outputdir/time_vs_clust$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsclust "\# cmax cmin c2sizeav c2sizestd nc2 csizeav csizestd nc clenav clenstd nc"
		}

	    }
	    "density_profile_calc" {

		# Initialize modes2d
		if { $profilenogrid } {
		    if { [catch { analyze bilayer_density_profile $hrange $nbins $beadtypes setgrid $mgrid $mgrid 0 setstray $stray_cut_off nogrid } ] } {
			mmsg::err $this "could not initialize density_profile_calc"
		    }
		} else {
		    if { [catch { analyze bilayer_density_profile $hrange $nbins $beadtypes setgrid $mgrid $mgrid 0 setstray $stray_cut_off } ] } {
			mmsg::err $this "could not initialize density_profile_calc"
		    }
		}

		#Initialize av_densities
		set thisbinlist 0.0
		unset thisbinlist
		for { set bn 0 } { $bn < $nbins } { incr bn } {
		    for { set bt 0 } { $bt < [expr 2*[llength $beadtypes]] } { incr bt } {
			lappend thisbinlist 0.0
		    }
		    lappend av_densities $thisbinlist
		    unset thisbinlist

		}
	    }
	    "localheights_calc" {
		set av_localheights 0
		unset av_localheights
		#Initialize av_localheights
		for { set bn 0 } { $bn < $localheightsnbins } { incr bn } {
		    lappend av_localheights 0
		}
	    }
	    "localorients_calc" {
		set av_localorients 0
		unset av_localorients
		#Initialize av_localorients
		for { set bn 0 } { $bn < $localorientsnbins } { incr bn } {
		    lappend av_localorients 0
		}
	    }
	    "pik1_calc" {
		for { set j 0 } { $j < $n_particles } { incr j } {
		    lappend all_particles $j
		}
		mmsg::debug $this "opening $outputdir/time_vs_pik1$suffix "

		if { [file exists "$outputdir/time_vs_pik1$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}
		set f_tvspik1 [open "$outputdir/time_vs_pik1$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvspik1 "\# Components of the total pressure tensor in row major order"
		    puts $f_tvspik1 "\# pxx pxy pxz pyx pyy pyz pzx pzy pzz"
		}
	    }
	    "pressure_calc" {
		mmsg::debug $this "Opening $outputdir/time_vs_pressure$suffix "

		if { [file exists "$outputdir/time_vs_pressure$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}

		set f_tvsp [open "$outputdir/time_vs_pressure$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsp "\# Components of the total pressure"
		    puts $f_tvsp "\# Time p_inst(NPT only) total ideal fene harmonic nonbonded "
		}
	    }
	    "box_len_calc" {
		mmsg::debug $this "opening $outputdir/time_vs_boxl$suffix "

		if { [file exists "$outputdir/time_vs_boxl$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}
		set f_tvsbl [open "$outputdir/time_vs_boxl$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsbl "\# Time boxx boxy boxz"
		}

	    }
	    "fluctuation_calc" {
		for { set r 0 } { $r < [expr $mgrid*($mgrid)] } { incr r } {
		    lappend av_pow 0.0
		}
		set rawmodes "[analyze modes2d setgrid $mgrid $mgrid 0 setstray $stray_cut_off]"
	    }
	    "flipflop_calc" {
		mmsg::debug $this "opening $outputdir/time_vs_flip$suffix "

		if { [file exists "$outputdir/time_vs_flip$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}

		set f_tvsflip [open "$outputdir/time_vs_flip$suffix" $iotype]
		set f_tvsflip [open "$outputdir/time_vs_flip$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsflip "\# Note: F(t) = N_up(t) + N_d(t)/N. See long membrane paper for details"
		    puts $f_tvsflip "\# Time flip_parameter"
		}
		set l_orients_start [ analyze get_lipid_orients setgrid $mgrid $mgrid 0 setstray $stray_cut_off ]
	    }
	    "orient_order_calc" {

		if { [file exists "$outputdir/time_vs_oop$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}
		mmsg::debug $this "opening $outputdir/time_vs_oop$suffix "
		set f_tvsoop [open "$outputdir/time_vs_oop$suffix" $iotype ]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsoop "\# S = 0.5*<3*(a_i.n)^2 -1>i"
		    puts $f_tvsoop "\# where a_i is the orientation of the ith lipid and n is the average bilayer normal"
		    puts $f_tvsoop "\# Time S"
		}
		

	    }
	    "stray_lipids_calc" {
		if { [file exists "$outputdir/time_vs_stray$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}
		mmsg::debug $this "opening $outputdir/time_vs_stray$suffix "
		set f_tvsstray [open "$outputdir/time_vs_stray$suffix" $iotype ]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsstray "\# The number of stray lipids vs time"
		    puts $f_tvsstray "\# Time num_strays"
		}
		set tmp [ analyze get_lipid_orients setgrid $mgrid $mgrid 0 setstray $stray_cut_off ]
	    }
	    "energy_calc" {

		if { [file exists "$outputdir/time_vs_energy$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}
		mmsg::debug $this "opening $outputdir/time_vs_energy$suffix "
		set f_tvsen [open "$outputdir/time_vs_energy$suffix" $iotype ]
		# First work out the names of components
		set raw [analyze energy]

		if { $newfile || $iotype == "w"} {
		    puts $f_tvsen "\# Components of the energy "
		    puts -nonewline $f_tvsen "\# Time total_en fene_en kinetic_en harmonic_en nonbonded_en"
		}
		unset av_components_en
		for { set k 0 } { $k < [llength $raw ] } { incr k } {
		    set tmp [lindex $raw $k]
		    set ntmp [llength $tmp]	
		    if { [ regexp "nonbonded" $tmp ] } {
			puts -nonewline $f_tvsen " [lrange $tmp 0 end-1]"
			lappend av_components_en 0.0
		    }
		}
		puts $f_tvsen ""

	    }
	    "distance_calc" {
		mmsg::debug $this "opening $outputdir/time_vs_distance$suffix "
		
		if { [file exists "$outputdir/time_vs_distance$suffix"] } {
		    set newfile 0
		} else {
		    set newfile 1
		}
		set f_tvsdist [open "$outputdir/time_vs_distance$suffix" $iotype]
               if { $newfile || $iotype == "w"} {
                   puts $f_tvsdist "\# Time Distance"
               }
		
	    }	    
	    "default" {
		mmsg::warn $this "unknown analysis flag [lindex $switches $i 0] $known_flags" 
	    }
	}
	
    }
    mmsg::debug $this "done"
    flush stdout
}

# ::std_analysis::finish_analysis --
#
# Cleanup by unsetting variables closing files etc.  If this is called
# then it should be OK to call setup_analysis to create an entirely
# new analysis setup.
#
proc ::std_analysis::finish_analysis {  } { 
    mmsg::send $this "finishing analysis"
    variable all_particles
    variable rawmodes
    variable av_pow
    variable av_pow_i 
    variable l_orients_start
    variable known_flags

    #File Streams
    variable f_tvspik1
    variable f_tvsp
    variable f_tvsbl
    variable f_tvsflip
    variable f_tvsoop
    variable f_tvsstray
    variable f_tvsen
    variable f_tvsclust
    variable f_tvsdist

    # Averaging
    variable av_dist
    variable av_dist_i

    variable av_localheights
    variable av_localorients

    variable av_densities
    variable av_densities_i

    variable av_clust
    variable av_clust_i
    
    variable av_pik1
    variable av_pik1_i 
    
    variable av_pressure
    variable av_pressure_i
    
    variable av_boxl
    variable av_boxl_i 
    
    variable av_flip 
    variable av_flip_i 
    
    variable av_oop 
    variable av_oop_i 
    
    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_fene_en
    variable av_harm_en
    variable av_nb_en
    variable av_en_i
    
    variable av_sizehisto
    variable av_sizehisto_i
    
    variable switches
    
    
    puts $switches	
    for { set i 0 } { $i < [llength $switches ] } { incr i } {
	puts [lindex $switches $i 0]
	puts -nonewline "Closing down: "
	switch [lindex $switches $i 0] {
	    "density_profile_calc" {
		puts "density_profile_calc"
		unset av_densities
	    }
	    "localheights_calc" {
		puts "localheights_calc"
		unset av_localheights
	    }
	    "localorients_calc" {
		puts "localorients_calc"
		unset av_localorients
	    }
	    "cluster_calc" {
		puts "cluster_calc"
		close $f_tvsclust
		set av_sizehisto 0
		set av_sizehisto_i 0
		set av_clust 0
		set av_clust_i 0
	    }
	    "pik1_calc" {
		puts "pik1_calc"
		catch { unset all_particles }
		close $f_tvspik1
		set av_pik1 0
		set av_pik1_i 0
	    }
	    "pressure_calc" {
		puts "pressure_calc"
		set av_pressure 0
		set av_pressure_i 0
		catch { unset all_particles }
		close $f_tvsp
	    }
	    "box_len_calc" {
		puts "box_len_calc"
		set av_boxl 0
		set av_boxl_i 0
		close $f_tvsbl 
	    }
	    "fluctuation_calc" {
		puts "fluctuation_calc"
		unset av_pow
		set av_pow_i 0
	    }
	    "flipflop_calc" {
		puts "flipflop_calc"
		set av_flip 0
		set av_flip_i 0
		close $f_tvsflip
	    }
	    "orient_order_calc" {
		puts "orient_order_calc"
		set av_oop 0
		set av_oop_i 0
		close $f_tvsoop
	    }
	    "stray_lipids_calc" {
		puts "stray_calc "
		close $f_tvsstray
	    }
	    "energy_calc" {
		puts "energy_calc"
		set av_components_en 0
		set av_total_en 0
		set av_kin_en 0
		set av_fene_en 0
		set av_harm_en 0
		set av_nb_en 0
		set av_en_i 0
		close $f_tvsen 		    
	    }
	    "default" {
		mmsg::warn $this "unknown analysis flag [lindex $switches $i 0] $known_flags"
	    }
	}
	
    }
    unset switches
    mmsg::debug $this "done"
}





