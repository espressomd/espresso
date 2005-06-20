

# ::std_analysis::analyze_clusters --
#
# Calculate the number of clusters.  Uses the "analyze aggregation"
# feature of espresso by Mehmet
#

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_clusters { printflag } {
    variable this
    mmsg::send $this "analyzing number of clusters"
    variable area_lipid
    variable av_clust
    variable av_clust_i
    variable av_sizehisto
    variable av_sizehisto_i

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
    set c2sizeav [::mathutils::average $c2sizedist]
    set c2sizestd [::mathutils::stdev $c2sizedist]
    set csizeav [::mathutils::average $sizedist]
    set csizestd [::mathutils::stdev $sizedist]
    set clenav [::mathutils::average $lendist]
    set clenstd [::mathutils::stdev $lendist]

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
    if { $printflag } {
	mmsg::send $this "i: $av_clust_i vals: $av_clust"
	mmsg::send $this  "cluster size av: [::mathutils::average $sizedist] std: [::mathutils::stdev $sizedist]"
	mmsg::send $this  "interface length av: [::mathutils::average $lendist] std: [::mathutils::stdev $lendist]"
	mmsg::send $this  "clusters: max: $cmax min: $cmin num: $cnum" 
    }
    set sizedist 0.0 
    set lendist 0.0
    set c2sizedist 0.0
    unset c2sizedist
    unset sizedist
    unset lendist

    mmsg::debug $this "done"

}
