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



proc calcObsAv { fileN ind } {
    # derives time-averages of the observables at the columns $ind in file $fileN,
    # returning '<amount of samples> { names (taken from header) } { averaged values }'
    set f [open $fileN "r"]
    gets $f tmp_line
    foreach i $ind { 
	set tmp "[lindex $tmp_line $i]"; lappend var $tmp
	set tmp "[lindex $tmp_line $i]_av"; lappend var_av $tmp
	eval set $tmp 0
    }
    while { [eof $f]==0 } { if { [gets $f tmp_line] > 0 } {
	    foreach i $ind j $var_av k $var {
		set tmp [lindex $tmp_line $i]; set $j [eval expr $$j + $tmp]; lappend $k $tmp
	    }
    } }
    close $f
    set N_av [llength $Temp]
    set res1 ""; set res2 ""; foreach i $var { lappend res1 $i }
    foreach j $var_av {	eval set $j [eval expr $$j / $N_av]; eval lappend res2 $$j }
    set res "$N_av \{ $res1 \} \{ $res2 \}"
    return $res

    puts -nonewline "$N_av samples => "
    foreach j $var_av {
	eval set $j [eval expr $$j / $N_av]; eval puts -nonewline \"$j = $$j \"
    }
    puts " "
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
    if { [llength $val]!=[llength $res] } { puts "\nWARNING: Not all values ($val) have been found in $what!\n" }
    return $res
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

proc plotObs { destinations } {
    set f [open "plotObsTmp.p" "w"]
    foreach destination $destinations {
	# Creates a gnuplot of the observables stored in $destination
	# assuming they are stored in the order 'time mindist re rg rh g123 Temp'
	# (as after a call to '[analysis ... [analyze energy kin]/$n_part/1.5]')
	if {[string first "KKG" $destination] == -1} {
	    puts $f "set xlabel \"time (LJ-time)\"; set ylabel \"$destination.a\""
	    puts $f "set nologscale xy"
	    puts -nonewline $f "plot \"$destination\" using 1:2 title \"mindist\" , "
	    puts -nonewline $f "\"$destination\" using 1:3 title \"R_e\" , "
	    puts -nonewline $f "\"$destination\" using 1:4 title \"R_g\" , "
	    puts -nonewline $f "\"$destination\" using 1:5 title \"R_h\" , "
	    puts -nonewline $f "\"$destination\" using 1:9 title \"Temp\"  \n"
	    puts $f "set out \"$destination.a.ps\""
	    puts $f "set terminal postscript color; replot"
	    puts $f "\# !lpr -Pthps18 \"$destination.a.ps\""
	    puts $f "set terminal x11"
	    puts $f " "
	    puts $f " "
	    puts $f "set xlabel \"time (LJ-time)\"; set ylabel \"$destination.b\""
	    puts $f "set logscale xy"
	    puts -nonewline $f "plot \"< tr '\{' ' ' < $destination\" using 1:6 title \"g1\" , "
	    puts -nonewline $f " \"< tr '\{' ' ' < $destination\" using 1:7 title \"g2\" , "
	    puts -nonewline $f " \"< tr '\{' ' ' < $destination\" using 1:8 title \"g3\"   \n"
	    puts $f "set out \"$destination.b.ps\""
	    puts $f "set terminal postscript color; replot"
	    puts $f "\# !lpr -Pthps18 \"$destination.b.ps\""
	    puts $f "set terminal x11"
	    puts $f " "
	    puts $f " "
	} else {
	    # Same as above, but for the (hopefully time-averaged) <re> and <rg> compared to
	    # the original KKG-values
	    puts $f "set xlabel \"monomers per chain\"; set ylabel \"$destination\""
	    puts $f "set logscale xy"
	    puts -nonewline $f "plot \"$destination\" using 3:6 title \"R_e\" w linespoints , "
	    puts -nonewline $f " \"$destination\" using 3:8 title \"R_e_KG\" w linespoints , "
	    puts -nonewline $f " \"$destination\" using 3:9 title \"R_g\" w linespoints , "
	    puts -nonewline $f " \"$destination\" using 3:11 title \"R_g_KG\" w linespoints  \n"
	    puts $f "set out \"$destination.ps\""
	    puts $f "set terminal postscript color; replot"
	    puts $f "\# !lpr -Pthps18 \"$destination.ps\""
	    puts $f "set terminal x11"
	    puts $f " "
	    puts $f " "
	    puts $f "set xlabel \"monomers per chain\"; set ylabel \"deviation in $destination\" "
	    puts $f "set logscale x; set nologscale y "
	    puts -nonewline $f "plot \"$destination\" using 3:7 title \"<R_e>-Re_KG\" , "
	    puts -nonewline $f " \"$destination\" using 3:10 title \"<R_g>-Rg_KG\" \n"
	    puts $f "set out \"$destination.div.ps\" "
	    puts $f "set terminal postscript color; replot "
	    puts $f "\# !lpr -Pthps18 \"$destination.div.ps\" "
	    puts $f "set terminal x11"
	    puts $f " "
	    puts $f " "
	}
    }
    close $f
    eval exec gnuplot plotObsTmp.p
    eval exec rm -f plotObsTmp.p
}

proc plotJoin { destinations final } {
    foreach destination $destinations {
	if {[string first "KKG" $destination] == -1} {
	    eval exec gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$destination.j.ps $destination.a.ps $destination.b.ps -c quit
	} else { 
	    eval exec gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$destination.j.ps $destination.ps $destination.div.ps -c quit    
	}
	lappend destinationsJ "$destination.j.ps"
    }
    eval exec gs -dNOPAUSE -sDEVICE=pswrite -sOutputFile=$final.1.ps $destinationsJ -c quit
    eval exec pstops '2:0L@.7(21cm,0)+1L@.7(21cm,14.85cm)' $final.1.ps $final.2.ps
}
