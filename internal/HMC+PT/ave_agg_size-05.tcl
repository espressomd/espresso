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

set filsup [open "super.txt" "w"]

set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]
set index 0

set startlist { 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 3000  3000 3000 3000 3000 3000 3000 3000 4000 4000 4000 4000 4000 5500 5500 5500 5500 5500  }

for {set i 0} { $i < $ensemble_num } {incr i} {
    set bjerrum [lindex $bjerrum_list $i]
    #create an array to store the number of aggregates
    set start [lindex $startlist $i]

    for {set k 0} {$k <= $n_poly} {incr k} {
	lappend num_agg 0
	lappend num_agg_tri 0
	lappend num_agg_mon 0
    }

    set counter 0
    set file [open "lb_[format %06.3f $bjerrum]/cionaggregation.txt" "r"]
    set filemax [open "lb_[format %06.3f $bjerrum]/maxaggregation.txt" "w"]
    puts "BJERRUM $bjerrum START $start"
    for { set config 0} { $config <= $end } { incr config $step} {
	set inp [gets $file]
	if { $config >= $start } {
	    incr counter
	    #puts "\n $config inp $inp"
	    set countmax 0
	    unset -nocomplain realaggs
	    for {set j 0 } { $j < [lindex $inp 11]} {incr j} {
		set cur_list [lindex $inp [expr 13 + $j]]
		set curr [llength $cur_list]
		set count 0
		set count_tri 0
		set count_mon 0
		for { set k 0} {$k < $curr } {incr k} {
		    if {[lindex $cur_list $k] < $n_poly } {
			incr count
		    } elseif { [lindex $cur_list $k] < [expr $n_poly + $num_trivalent ] } {
			incr count_tri
		    } else {
			incr count_mon
		    }
		}
		
		for { set k [expr $curr-1]} {$k > -1 } {incr k -1} {
		    if { [lindex $cur_list $k] >= $n_poly } { set cur_list [lreplace $cur_list $k $k ]}
		}
		if {[llength $cur_list] > 0 } {lappend realaggs $cur_list}

		if {$count > 0} {
		    if {$countmax < $count } {set countmax $count}
		    lset num_agg $count [expr [lindex $num_agg $count] +1] 
		}
		if {$count_tri > 0} {
		    lset num_agg_tri $count [expr [lindex $num_agg_tri $count] + $count_tri] 
		}
		if {$count_mon > 0} {
		    lset num_agg_mon $count [expr [lindex $num_agg_mon $count] + $count_mon ] 
		}

	    }

	    set countmin $n_poly
	    set countavg 0
	    set countstd 0
	    set aggnum [llength $realaggs]
	    for {set jj 0} { $jj < $aggnum} {incr jj} {
		set clen  [llength [lindex $realaggs $jj] ]
		if { $clen < $countmin } {set countmin $clen}
		set countavg [expr $countavg + $clen]
		set countstd [expr $countstd + $clen * $clen]
	    }
	    set countavg [expr double($countavg) / $aggnum]
	    set countstd [expr sqrt(double($countstd) / $aggnum - $countavg * $countavg) ]

	    puts $filemax "Run $config MAX $countmax MIN $countmin AVG $countavg STD $countstd AGG_NUM $aggnum AGGREGATES $realaggs"  
#Run 0  MAX 1 MIN 1 AVG 1.000000 STD 0.000244 AGG_NUM 61 AGGREGATES"

	    
	    #puts "$inp"
	}
    }
    #puts "NUM_AGG $num_agg"

    set n_avg 0
    set m_avg 0
    set stan 0
    set tot 0
    for {set k 0} {$k <= $n_poly} {incr k} {
	set n_avg [expr $n_avg + [lindex $num_agg $k] * ($k)]
	set m_avg [expr $m_avg + [lindex $num_agg $k] * ($k) * ( $k)]
	set tot [expr [lindex $num_agg $k] + $tot ]
    }

    lset num_agg_tri 0 [expr  double([lindex $num_agg_tri 0]) / $counter / $num_trivalent]
    if { $num_monovalent >0} {lset num_agg_mon 0 [expr  double([lindex $num_agg_mon 0]) / $counter / $num_monovalent]}
    for {set k 1} {$k <= $n_poly} {incr k} {
	if {[lindex $num_agg $k] != 0} {
	    lset num_agg_tri $k [expr 1. - double([lindex $num_agg_tri $k ]) / [lindex $num_agg $k] / $p_length * $c_ion_val / $k]
	    if {  $monovalentcon > 0} {lset num_agg_mon $k [expr 1. - double([lindex $num_agg_mon $k ]) / [lindex $num_agg $k] / $p_length * $c_ion_val / $k]}
	    lset num_agg $k [expr double([lindex $num_agg $k] ) / $counter * $k / $n_poly]
	}
    }
    set wa [expr double($m_avg)/$n_avg]
    set na [expr double($n_avg) / $tot]
    set pdi [expr $wa / $na]
    #puts $filsup "BJERRUM\t[format %5.2f $bjerrum]\tWA\t$wa\tNA $na PDI $pdi DIST $num_agg TRI $num_agg_tri MON $num_agg_mon"
    puts -nonewline $filsup "BJERRUM\t[format %05.3f $bjerrum]\tWA\t[format %010.5f $wa]\tNA\t[format %010.5f $na]\tPDI\t[format %010.5f $pdi]\tDIST\t"
    for {set ij 0} { $ij < [llength $num_agg]} {incr ij} {puts -nonewline $filsup "[lindex $num_agg $ij]\t" }
    puts -nonewline $filsup "TRI\t"
    for {set ij 0} { $ij < [llength $num_agg_tri]} {incr ij} {puts -nonewline  $filsup "[lindex $num_agg_tri $ij]\t" }
    puts -nonewline $filsup "MON\t"
    for {set ij 0} { $ij < [llength $num_agg_mon]} {incr ij} {puts -nonewline $filsup "[lindex $num_agg_mon $ij]\t" }
    puts $filsup ""

    set num_agg_$i $num_agg
    #puts "HE [expr \$num_agg_$i ] "
    
    unset num_agg 
    unset num_agg_tri
    unset num_agg_mon
    close $file
    close $filemax
}


set file [open "distaggregate.txt" "w"]
for {set k 1} {$k <= $n_poly} {incr k} {
    puts -nonewline $file "Agg [expr $k] "
    for {set i 0} { $i < $ensemble_num } {incr i} {
	puts -nonewline $file " [lindex [expr \$num_agg_$i ] $k]"
    }
    puts -nonewline $file "\n"
}
puts $file "\n"


exit 0
