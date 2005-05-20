#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#  
#############################################################
#                                                           #
# statistics.tcl                                            #
# ==============                                            #
#                                                           #
# Some scripts for statistical analysis.                    #
#                                                           #
# Created:       01.04.2003 by BAM                          #
#                                                           #
#############################################################



proc average_checkpoint { origin what } {
    # reads in a checkpoint, executing analysis-command 'what' after each file, averaging it in the end
    if { [file exists "$origin.chk"] } { set chk [open "$origin.chk" "r"] 
    } elseif { [file exists "$origin"] } { set chk [open "$origin" "r"] 
    } else { puts "ERROR: Could not find checkpoint-list $origin!\nAborting..."; exit }
    while { [eof $chk]==0 } { if { [gets $chk source] > 0 } {
	if { [string compare [lindex [split $source "."] end] "gz"]==0 } { set f [open "|gzip -cd $source" r]
	} else { set f [open "$source" "r"] }
	while { [blockfile $f read auto] != "eof" } {}
	puts -nonewline "."; flush stdout; # puts "read $source"
	close $f
    } }
    close $chk
}



proc calcObsAv { fileN ind { startJ 0 } } {
    # derives time-averages of the observables at the columns $ind in file $fileN,
    # returning '<amount of samples> { names (taken from header) } { averaged values } { errors }',
    # skip the first <startJ>-lines before starting to average
    set f [open $fileN "r"]
    gets $f tmp_line
    foreach i $ind { 
	set tmp "[string map { - _ } [lindex $tmp_line $i]]";     lappend var $tmp
	set tmp "[string map { - _ } [lindex $tmp_line $i]]_av";  lappend var_av $tmp; eval set $tmp 0
	set tmp "[string map { - _ } [lindex $tmp_line $i]]_av2"; lappend var_av2 $tmp; eval set $tmp 0
	
    }
    set N_av 0; set imax [lindex [lsort $ind] end]
    # skip first $startJ-lines before starting to average
    for {set i 0} {$i < $startJ} {incr i} { gets $f tmp_line }
    # now average the remainders of the file
    while { [eof $f]==0 } { if { [gets $f tmp_line] > 0 } {
	if {[llength $tmp_line] <= $imax } { puts "File corrupted, current line too short (got: '$tmp_line', but need at least $imax+1 entries)!" 
	} else {
	    foreach i $ind  j $var_av  j2 $var_av2  k $var { 
		set tmp [lindex $tmp_line $i]; set $j [eval expr $$j + $tmp]; set $j2 [eval expr $$j2 + $tmp*$tmp] }
	    incr N_av
	}
    } }
    close $f
    set res1 ""; set res2 ""; set res3 ""; foreach i $var { lappend res1 $i }
    if { $N_av == 1 } {
	foreach j $var_av {	eval set $j [eval expr $$j / $N_av]; eval lappend res2 $$j }
    } elseif { $N_av > 1 } {
	foreach k $var j $var_av  j2 $var_av2 {
	    eval set $j [eval expr $$j / $N_av]; eval lappend res2 $$j
	    eval set $j2 [eval expr sqrt(abs($$j2/$N_av - $$j*$$j))]; eval lappend res3 $$j2 
	}
    }
    set res "$N_av \{ $res1 \} \{ $res2 \} \{ $res3 \}"
    return $res
}

proc findObsAv { val what } {
    # extracts the values in the tcl-list $val from $what (the output-list from 'calcObsAv')
    # returning them as a tiny tcl-list
    lappend res [lindex $what 0]; set tmp1 [lindex $what 1]; set tmp2 [lindex $what 2]
    foreach j $val {
	for { set i 0 } { $i < [llength $tmp1] } { incr i } {
	    if { [lindex $tmp1 $i] == $j } { lappend res [lindex $tmp2 $i] }
	}
    }
    if { [llength $val]!=[llength $res]-1 } { puts "\nWARNING: Not all values ($val) have been found in $what!\n" }
    return $res
}

proc calcObAv { fileN ind { startJ 0 } } {
    # does the same as 'calcObsAv', but for one observable only, hence returning only its averaged value
    if { [llength $ind]!=1 } { puts "\nWARNING: Parameter '$ind' is too long - use 'calcObsAv' to average multiple observables!"; exit }
    set what [calcObsAv $fileN $ind $startJ]
    return [lindex $what 2]
}

proc calcObErr { fileN ind { startJ 0 } } {
    # same as 'calcObAv', but returns the error of that value
    if { [llength $ind]!=1 } { puts "\nWARNING: Parameter '$ind' is too long - use 'calcObsAv' to average multiple observables!"; exit }
    set what [calcObsAv $fileN $ind $startJ]
    return [lindex $what 3]
}

proc nameObsAv { fileN names { startJ 0 } } {
    # does the same as 'calcObsAv', but expects the observables' column-names rather than their column-positions
    set f [open $fileN "r"]; set ind1 ""; set ind2 ""
    gets $f tmp_line
    for { set j 0 } { $j<[llength $names] } { incr j } {
	for { set i 0 } { $i<[llength $tmp_line] } { incr i } {
	    if { "[lindex $tmp_line $i]" == "[lindex $names $j]" } { lappend ind1 $i; lappend ind2 $j }
	}
    }
    close $f
    if { $ind2=="" || $ind1=="" } { puts "\nERROR: None of the observables were found (you were looking for '$names')!"; exit }
    if { [llength $names]!=[llength $ind2] } { puts "\nWARNING: Only [llength $ind2] of [llength $names] parameters have been found in $fileN!"; exit }
    set what [calcObsAv $fileN $ind1 $startJ]
    return [concat [lindex $what 0] [lindex $what 2] [lindex $what 3]]
}


proc replObsAv { fileN what } {
    # replaces the values for 're' and 'rg' in $fileN by their time-averages
    # using the output-list from 'calcObsAv', writing the result to '$fileN.av'
    set reD -1; set rgD -1; set reO -1; set rgO -1
    set fin  [open $fileN "r"]
    set fout [open "$fileN.av" "w"]
    gets $fin tmp_line
    for { set i 0 } { $i < [llength $tmp_line] } { incr i } {
	if {[lindex $tmp_line $i]=="re"} { set reI $i; puts -nonewline $fout "re_av " 
	} elseif {[lindex $tmp_line $i]=="rg"} { set rgI $i; puts -nonewline $fout "rg_av "
	} elseif {[lindex $tmp_line $i]=="re-reKG"} { set reD $i; puts -nonewline $fout "re_av-reKG reKG"
	} elseif {[lindex $tmp_line $i]=="rg-rgKG"} { set rgD $i; puts -nonewline $fout "rg_av-rgKG rgKG"
	} elseif {[lindex $tmp_line $i]=="reKG"} { set reO $i 
	} elseif {[lindex $tmp_line $i]=="rgKG"} { set rgO $i 
	} else { puts -nonewline $fout "[lindex $tmp_line $i] " }
    }
    puts $fout " "
    set tmp_what [lindex [lindex $what 0] 1]
    for { set i 0 } { $i < [llength $tmp_what] } { incr i } {
	if {[lindex $tmp_what $i]=="re"} { set reW $i } elseif {[lindex $tmp_what $i]=="rg"} { set rgW $i }
    }
    set posW 0
    while { [eof $fin]==0 } { if { [gets $fin tmp_line] > 0 } {
	set tmp_what [lindex [lindex $what $posW] 2]
	for { set i 0 } { $i < [llength $tmp_line] } { incr i } {
	    if { $i==$reI } { puts -nonewline $fout "[lindex $tmp_what $reW] "
	    } elseif { $i==$rgI } { puts -nonewline $fout "[lindex $tmp_what $rgW] "
	    } elseif { $i==$reD } { 
		puts -nonewline $fout "[expr ([lindex $tmp_what $reW]-[lindex $tmp_line $reO])/[lindex $tmp_line $reO]] "
		puts -nonewline $fout "[lindex $tmp_line $reO] "
	    } elseif { $i==$rgD } { 
		puts -nonewline $fout "[expr ([lindex $tmp_what $rgW]-[lindex $tmp_line $rgO])/[lindex $tmp_line $rgO]] "
		puts -nonewline $fout "[lindex $tmp_line $rgO] "
	    } elseif { $i!=$reO && $i!=$rgO } { puts -nonewline $fout "[lindex $tmp_line $i] " }
	}
	puts $fout " "; incr posW
    } }
    close $fin; close $fout
    return $posW
}




#
# gnuplot
# -------
# 
# Some auxiliary functions for creating gnuplot-outputs
# of the simulation data.
#
#############################################################

proc plotObs { destinations what {p1 NA} {p2 NA} {p3 NA} {p4 NA} {p5 NA} {p6 NA} {p7 NA} {p8 NA} {p9 NA} {p10 NA} } {
    # Creates a gnuplot of the data in $destination at positions $what.
    # Syntax: 'plotObs <data> { x:y1 ... x:yn } [titles { "title.y1" ... "title.yn" }] [labels { "xlabel" "ylabel" }] [scale <gnuplot-scale>] [out <out>] [cmd <gnuplot-command>]'
    set param [list $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10]
    for {set i 0} {$i < [llength $what]} {incr i} { lappend titles "Data $i" }
    set labels [list "x-values" "$destinations"]; set scale "nologscale xy"; set out [lindex $destinations end]; set cmd ""
    for {set i 0} {$i < 10} {incr i} {
	switch [lindex $param $i] {
	    "titles" { incr i; set titles [lindex $param $i]
		if { [llength $titles] > [llength $what] } { set titles [lrange $titles 0 [llength $what]] }
		if { [llength $titles] < [llength $what] } { 
		    for {set j [llength $titles]} {$j < [llength $what]} {incr j} {lappend titles "Data $j"}
		}  }
	    "labels" { incr i; set labels [lindex $param $i]
		if { [llength $labels]!=2 } { set labels [list "[lindex $labels 0]" "$destinations"] }  }
	    "scale"  {incr i; set scale [lindex $param $i]  }
	    "out"    {incr i; set out [lindex $param $i]    }
	    "cmd"    {incr i; set cmd [lindex $param $i]    }
	    default { if { [lindex $param $i]!="NA" } {
		puts "The parameter set you supplied ($param) does not seem to be valid (stuck at: [lindex $param $i])!\nAborting...\n"; exit }
	    }
	}
    }
    if {[llength $destinations] < [llength $what]} {
	for {set i [llength $destinations]} {$i < [llength $what]} {incr i} { lappend destinations [lindex $destinations end] }
    }
    set f [open "plotObsTmp.p" "w"]
    puts $f "set out \"$out.ps\""
    puts $f "set terminal postscript color"
    puts $f "set xlabel \"[lindex $labels 0]\"; set ylabel \"[lindex $labels 1]\""
    puts $f "set $scale"
    puts -nonewline $f "plot "
    for {set i 0} {$i < [expr [llength $what]-1]} {incr i} {
	puts -nonewline $f "\"[lindex $destinations $i]\" using [lindex $what $i] title \"[lindex $titles $i]\" , "
    }
    puts -nonewline $f "\"[lindex $destinations end]\" using [lindex $what $i] title \"[lindex $titles $i]\"  \n"
    if { $cmd != "" } { puts $f $cmd }
    # puts $f "replot"
    puts $f "\# !lpr -Pthps18 \"$out.ps\""
    puts $f "set terminal x11"
    close $f
    catch { eval exec gnuplot plotObsTmp.p }
    # eval exec gnuplot plotObsTmp.p
    eval exec rm -f plotObsTmp.p
}

proc plotJoin { sources final } {
    eval exec  gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$final.A.ps $sources -c quit
    catch { eval exec pstops "2:0L@.7(21cm,0)+1L@.7(21cm,14.85cm)" $final.A.ps $final.B.ps }
}
