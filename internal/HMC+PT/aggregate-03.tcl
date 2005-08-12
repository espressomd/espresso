#!/bin/sh
source ./hmc-initialize.tcl
source ./hmc-setup.tcl

set f [open "input.txt" "r"]
while { [eof $f] == 0 } {
    gets $f temp;
    if {[lindex $temp 0] == "monovalentcon" } { set monovalentcon [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "trivalentcon"} { set trivalentcon [lindex $temp 1]} 
}
close $f


#this is the total number of pe's + 3 valent c-ions
set num_trivalent [expr  int ($n_poly * $p_length / $c_ion_val *  $trivalentcon / 100 ) ]
set num_monovalent [expr int ($n_poly * $p_length / 1 * $monovalentcon / 100  )]

# -1 because this is actually not the total aggregate number but max_aggregate id
set total_agg [expr $n_poly + $num_trivalent + $num_monovalent -1 ]

set ljr_cut  10.0
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
inter 0 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 0 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 1 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 1 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 2 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0

set f [open "|gzip -cd equilibrized.gz" "r"]
while { [blockfile $f read auto] != "eof" } {}
close $f

# set the topology
for {set i 0} { $i < $n_poly} {incr i} {
    set topo_chain 0 
    for {set j 0} { $j < $p_length} {incr j} {
	lappend topo_chain [expr $i * $p_length + $j]
    }
    lappend topo $topo_chain
    #puts "$topo_chain \n $topo"
}

#for the trivalent counterions set molecule type to n_poly+1
for {set i 0} { $i < [expr $num_trivalent]} {incr i} {
    set topo_chain  0 
    lappend topo_chain [expr $n_mono + $i] 
    lappend topo $topo_chain
}

#for the monovalent counterions set molecule type to n_poly+1
for {set i 0} { $i < [expr $num_monovalent]} {incr i} {
    set topo_chain  0 
    lappend topo_chain [expr $n_mono + $num_trivalent + $i] 
    lappend topo $topo_chain
}


#puts "$topo"

eval analyze set $topo
analyze set "topo_part_sync"

set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]
set index 0

for {set index 0} { $index < $ensemble_num} {incr index} {
    set bjerrum [lindex $bjerrum_list $index]
    #if { $bjerrum <= 1.500} {continue}

    puts "*****************************\n BJERRUM $bjerrum"
    set f [open "lb_[format %06.3f $bjerrum]/cionaggregation.txt" a]
    #########################################################################
    for { set config $start } { $config <= $end } { incr config $step} {
	#exec less lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz > coordinates
	#set file [open "coordinates" "r"]
	#set file [open "|less  lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
	#set file [open "|bzip2 -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].bz2" "r"]
	set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
	#set file [open "lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config]" "r"]
	while { [blockfile $file read auto] != "eof" } {}
#	close $file
	if {[catch {close $file} err]} {
	    if {![regexp "decompression OK, trailing garbage ignored" $err]} {
		error $err
	    }
	}

	set result [analyze charge_aggregation [expr $bjerrum] 0 $total_agg]
	#set result [analyze aggregation 3. 0 $total_agg]
	if {[lindex $argv 4] == "file" } {
	    puts $f "Run $config $result" 
	} else {
	    puts "$result"
	}
    }
    close $f
}


exit 0