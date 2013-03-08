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


# ::mbtools::analysis::analyze_clusters --
#
# Calculate the number of clusters.  Uses the "analyze aggregation"
# feature of espresso by Mehmet
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::clusters {
    variable area_lipid
    variable av_clust { 0 0 0 0 0 0 0 0 0 0 0 }
    variable av_clust_i 0

    variable av_sizehisto 0
    variable av_sizehisto_i 0

    variable f_tvsclust

    variable verbose

    namespace export setup_clusters
    namespace export analyze_clusters
    namespace export printav_clusters
    namespace export resetav_clusters

}



proc ::mbtools::analysis::clusters::resetav_clusters { } {
    variable av_clust 
    variable av_clust_i 
    variable av_sizehisto
    variable av_sizehisto_i

    set av_clust {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_clust_i 0
    set av_sizehisto 0
    set av_sizehisto_i 0
}

proc ::mbtools::analysis::clusters::printav_clusters { } {
    variable av_clust_i
    variable av_clust
    variable av_sizehisto
    variable av_sizehisto_i
    variable f_tvsclust
    global ::mbtools::analysis::time
    global ::mbtools::analysis::outputdir

    if { $av_clust_i > 0 } {
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
	::mmsg::warn [namespace current] "can't print average clusters"
	flush stdout
    }
    flush $f_tvsclust
}


proc ::mbtools::analysis::clusters::setup_clusters { args } {
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    global ::outputdir
    variable verbose
    variable area_lipid
    variable f_tvsclust
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_clust$suffix "		    
    
    set options {
	{alipid.arg  1.29 "area per lipid" }
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_fluctuations gridm:straycutoff "
    array set params [::cmdline::getoptions args $options $usage]

    set area_lipid $params(alipid)
    set verbose $params(verbose)

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


proc ::mbtools::analysis::clusters::analyze_clusters { } {
    variable area_lipid
    variable av_clust
    variable av_clust_i
    variable av_sizehisto
    variable av_sizehisto_i
    variable verbose

    ::mmsg::send [namespace current] "analyzing number of clusters"

    set topo [analyze set]
    set nmols [llength $topo]
    set ntypes [lindex [lindex $topo [expr $nmols -1]] 0]

    set k 0
    for { set i 0 } { $i < [llength $topo] } { incr i } {
	set typetmp [lindex [lindex $topo $i] 0]
	if { $typetmp != $k } {
	    lappend composition $i
	    incr k
	}
    }
    lappend composition $nmols
    puts $composition
    set result [analyze aggregation 1.1 0 [expr [lindex $composition 0] -1]]
    
    set cmax [lindex $result 1]
    set cmin [lindex $result 3]
    set cnum 0
    for { set i 11 } { $i < [llength $result] } { incr i } {
	set csize [llength [lindex $result $i]]
	if { $csize > 1 } {
	    incr cnum
	    lappend c2sizedist $csize
	}
	
	set carea [expr $csize*$area_lipid]
	set cradius [expr sqrt($carea/(1.0*3.141592654))]
	set clength [expr 2.0*$cradius*3.141592654 ]
	lappend sizedist $csize
	lappend lendist $clength
    }
    set c2sizeav [::mbtools::utils::average $c2sizedist]
    set c2sizestd [::mbtools::utils::stdev $c2sizedist]
    set csizeav [::mbtools::utils::average $sizedist]
    set csizestd [::mbtools::utils::stdev $sizedist]
    set clenav [::mbtools::utils::average $lendist]
    set clenstd [::mbtools::utils::stdev $lendist]

    # The maximum cluster size
    lset av_clust  0 [expr [lindex $av_clust 0] + $cmax]
    # The minimum cluster size
    lset av_clust  1 [expr [lindex $av_clust 1] + $cmin]
    # Average size of clusters including those of size 2 or greater
    lset av_clust  2 [expr [lindex $av_clust 2] + $c2sizeav]
    # Std deviation of clusters including those of size 2 or greater
    lset av_clust  3 [expr [lindex $av_clust 3] + $c2sizestd]
    # Number of clusters of size 2 or greater
    lset av_clust  4 [expr [lindex $av_clust 4] + [llength $c2sizedist]]
    # Total average cluster size
    lset av_clust  5 [expr [lindex $av_clust 5] + $csizeav]
    # Total cluster size std deviation
    lset av_clust  6 [expr [lindex $av_clust 6] + $csizestd]
    # Total number of clusters
    lset av_clust  7 [expr [lindex $av_clust 7] + [llength $sizedist]]
    # Length of the interface
    lset av_clust  8 [expr [lindex $av_clust 8] + $clenav]
    # std of the interface length
    lset av_clust  9 [expr [lindex $av_clust 9] + $clenstd]
    # Number of clusters for which we calculated the length
    lset av_clust  10 [expr [lindex $av_clust 10] + [llength $lendist]]

    #	puts $result
    
    #	puts "[llength $av_sizehisto] , $cmax"
    for { set i [llength $av_sizehisto] } { $i <= $cmax } { incr i } {
	lappend av_sizehisto 0 
    }
    
    #	puts "[llength $av_sizehisto] , $cmax, [llength $sizedist]"
    for { set i 0 } { $i < [llength $sizedist] } { incr i } {
	set bin [expr [lindex $sizedist $i] -1]
	lset av_sizehisto $bin [expr [lindex $av_sizehisto $bin ] + 1]	    
	#	    puts "$bin [lindex $av_sizehisto $bin] "
    }
    incr av_sizehisto_i
    
    
    
    incr av_clust_i
    if { $verbose } {
	::mmsg::send [namespace current] "i: $av_clust_i vals: $av_clust"
	::mmsg::send [namespace current]  "cluster size av: [::mbtools::utils::average $sizedist] std: [::mbtools::utils::stdev $sizedist]"
	::mmsg::send [namespace current]  "interface length av: [::mbtools::utils::average $lendist] std: [::mbtools::utils::stdev $lendist]"
	::mmsg::send [namespace current]  "clusters: max: $cmax min: $cmin num: $cnum" 
    }
    set sizedist 0.0 
    set lendist 0.0
    set c2sizedist 0.0
    unset c2sizedist
    unset sizedist
    unset lendist

    ::mmsg::debug [namespace current] "done"

}
