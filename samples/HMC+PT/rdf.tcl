#!/bin/sh
source ./hmc-initialize.tcl
source ./hmc-setup.tcl

set ljr_cut  10.0
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
inter 0 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 1 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0

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

#for the cions set molecule type to n_poly+1
set topo_chain  1 
for {set i 0} { $i < $n_ci } {incr i} {
    lappend topo_chain [expr $n_mono + $i] 
}
lappend topo $topo_chain
#puts "$topo\n"

eval analyze set $topo
analyze set "topo_part_sync"

set startlist { 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 3000  3000 3000 3000 3000 3000 3000 3000 4000 4000 4000 4000 4000 5500 5500 5500 5500 5500  }

set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]
set index 0
puts "I AM HERE  $lb_add $fac_lb $ensemble_num"
for {set index 0} { $index < $ensemble_num} {incr index} {
    set bjerrum [lindex $bjerrum_list $index]
    set start [lindex $startlist $index]

    if {(($bjerrum == 1.50) || ($bjerrum == 1.60) || ($bjerrum == 1.70) || ($bjerrum == 1.80) || ($bjerrum == 1.90) )} {} else {continue}
    #if {(($bjerrum < 1.940) || ($bjerrum == 1.80) || ($bjerrum >= 2.00))} {continue}
    set count 0
    unset -nocomplain totalrdf
    for {set i 0} {$i < 100 } {incr i} {
	lappend totalrdf {0. 0. 0. 0.}
    }

    
    puts "*****************************\n BJERRUM $bjerrum"
    #set f [open "lb_[format %06.3f $bjerrum]/aggregation.txt" a]
    #########################################################################
    for { set config $start } { $config <= $end } { incr config $step} {
	#puts  "Read checkpoint.[format %05d $config] - create pdb\r" 
	incr count 
 	puts  -nonewline "Read checkpoint.[format %05d $config]\r" 
	flush stdout
	set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
	while { [blockfile $file read auto] != "eof" } {}
	if {[catch {close $file} err]} {
	    if {![regexp "decompression OK, trailing garbage ignored" $err]} {
		error $err
	    }
	}

	set min_contact_num [expr $p_length/2]
	set contact_dis [expr $bjerrum * 3]
	set num_bins 100
	set rdf_mind 1.
	set rdf_maxd [expr $box_l/2]
	set s_molid 0
	set f_molid [expr $n_poly -1]
	#set result [analyze my_aggregation $contact_dis 0 [expr $n_poly-1] $min_contact_num ]
	#set rdf [analyze my_rdfchain $rdf_mind $rdf_maxd $num_bins $contact_dis $min_contact_num $s_molid $f_molid $s_molid $n_poly $p_length]
	set rdf [analyze rdfchain $rdf_mind $rdf_maxd $num_bins $s_molid $n_poly $p_length]
	#puts "$rdf"
	if {[lindex $argv 4] == "file" } {
	    #puts $f "Run $config $result" 
	} else {
	    #puts "$result"
	    for {set i 0} {$i < 100 } {incr i} {
		for {set j 0} { $j < 4} {incr j} {
		    lset totalrdf $i $j [expr [lindex $totalrdf $i $j ] + [lindex $rdf $i $j]]
		}
	    }
	}
    }
    #close $f

#NORMALIZE 
    for {set i 0} {$i < 100 } {incr i} {
	for {set j 0} { $j < 4} {incr j} {
	    lset totalrdf $i $j [expr [lindex $totalrdf $i $j ] / $count]
	}
    }

    set frdf [open "lb_[format %06.3f $bjerrum]/rdf.txt" "w"]
    for {set i 0} {$i < 100 } {incr i} {
	puts $frdf "[lindex $totalrdf $i]"
    }
    close $frdf
}


exit 0
