#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;

puts "----------------------------------------------"
puts "- Testcase analysis.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"
set errf [lindex $argv 1]

proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	set f [open $errf "w"]
	puts $f "not compiled in: $feature"
	close $f
	exit -42
    }
}

proc rewrite {in out} {
    global observables listables

    puts -nonewline "Importing checkpoints from '$in' and re-writing them to '$out'... "; flush stdout
    if { [file exists "$in.chk"] } { set chk [open "$in.chk" "r"] 
    } elseif { [file exists "$in"] } { set chk [open "$in" "r"] 
    } else { puts "ERROR: Could not find checkpoint-list $in!\nAborting..."; exit }
    set i 0; set sys_obs [list box_l gamma periodicity skin temperature time time_step]
    while { [eof $chk]==0 } { if { [gets $chk source] > 0 } {
	if { [string compare [lindex [split $source "."] end] "gz"]==0 } { set f [open "|gzip -cd $source" r]
	} else { set f [open "$source" "r"] }
	while { [blockfile $f read auto] != "eof" } {}
	puts -nonewline "."; flush stdout; close $f

	analyze set chains 0 20 30; set observables ""; set listables ""
	set mindist1 [analyze mindist]; set mindist2 [analyze mindist 0 0]; lappend observables $mindist1 $mindist2
	set nbhood [analyze nbhood 13 2.5]; lappend listables $nbhood
	set distto [analyze distto 13]; lappend observables $distto
	set energy [analyze energy total]; lappend observables $energy
	set pressure [analyze pressure total]; lappend observables $pressure
	set re [analyze re]; set re_av [analyze <re>]; lappend listables $re $re_av
	set rg [analyze rg]; set rg_av [analyze <rg>]; lappend listables $rg $rg_av
	set rh [analyze rh]; set rh_av [analyze <rh>]; lappend listables $rh $rh_av
	set idf [analyze internal_dist]; set idf_av [analyze <internal_dist>]; lappend listables $idf $idf_av
	set bdf [analyze bond_dist index 13]; set bdf_av [analyze <bond_dist> index 13]; lappend listables $bdf $bdf_av
	set bondl [analyze bond_l]; set bondl_av [analyze <bond_l>]; lappend listables $bondl $bondl_av
	set gff [analyze formfactor 1 10 10]; set gff_av [analyze <formfactor> 1 10 10]; lappend listables $gff $gff_av
	set g1v [analyze <g1>]; set g2v [analyze <g2>]; set g3v [analyze <g3>]; lappend listables $g1v $g2v $g3v
	checkpoint_set "$out.[eval format %02d $i].gz" 1 "observables listables" "all" "$sys_obs" "-"
	incr i
    } }
    close $chk; puts " Done."
}


if { [setmd n_nodes] == 3 || [setmd n_nodes] == 5 || [setmd n_nodes] == 6 || [setmd n_nodes] == 7} {
    puts "Testcase analysis.tcl does not run on 3,5,6 or 7 nodes"
    exec rm -f $errf
    exit 0
}

set epsilon  1e-4
setmd temp   0.0
setmd gamma  0.0
#set tcl_precision  6

### rewrite analysis_system.data analysis_system.data2; exit 0
### rewrite analysis_system_BB.data analysis_system.data; exit 0

if { [catch {
    puts -nonewline "Reading the checkpoint... "; flush stdout
    checkpoint_read analysis_system.data; puts " Done."

    analyze set chains 0 20 30
    set mindist1 "analyze mindist"; set mindist2 "analyze mindist 0 0"; lappend get_observables $mindist1 $mindist2
    set nbhood "lsort -integer \[analyze nbhood 13 2.5\]"; lappend get_listables $nbhood
    set distto "analyze distto 13"; lappend get_observables $distto
    set energy "analyze energy total"; lappend get_observables $energy
    set pressure "analyze pressure total"; lappend get_observables $pressure
    set re "analyze re"; set re_av "analyze <re>"; lappend get_listables $re $re_av
    set rg "analyze rg"; set rg_av "analyze <rg>"; lappend get_listables $rg $rg_av
    set rh "analyze rh"; set rh_av "analyze <rh>"; lappend get_listables $rh $rh_av
    set idf "analyze internal_dist"; set idf_av "analyze <internal_dist>"; lappend get_listables $idf $idf_av
    set bdf "analyze bond_dist index 13"; set bdf_av "analyze <bond_dist> index 13"; lappend get_listables $bdf $bdf_av
    set bondl "analyze bond_l"; set bondl_av "analyze <bond_l>"; lappend get_listables $bondl $bondl_av
    set gff "analyze formfactor 1 10 10"; set gff_av "analyze <formfactor> 1 10 10"; lappend get_listables $gff $gff_av
    set g1v "analyze <g1>"; set g2v "analyze <g2>"; set g3v "analyze <g3>"; lappend get_listables $g1v $g2v $g3v

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    
    # to ensure force recalculation
    invalidate_system
    integrate 0

    foreach obs $observables get_obs $get_observables {
	set rel_error [expr abs(([eval $get_obs] - $obs)/$obs)]
	puts "relative deviations upon evaluating '$get_obs': $rel_error  ([eval $get_obs] / $obs)"
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error too large upon evaluating '$get_obs'  ([eval $get_obs] / $obs)"
	}
    }

    foreach lst $listables get_lst $get_listables {
	set rel_max 0; set abs_max 0; set absflag 0; set maxi "-"; set maxj "-"
	foreach i $lst j [eval $get_lst] {
	    if { [string first "analyze formfactor" "$get_lst"]==0 || [string first "analyze <formfactor>" "$get_lst"]==0 } { 
		if { [expr [lindex $i 0]-[lindex $j 0]] < $epsilon } { 
		    set i [lindex $i 1]; set j [lindex $j 1] 
		} else { error "different x-coordinates upon comparing '$get_lst'" }
	    }
	    if { $i!=0 && $j!=0 } { set rel_error [expr abs(($j - $i)/$i)] } else { set rel_error -1; set absflag 1 }
	    set abs_error [expr abs($i-$j)]
	    if { $rel_error > $epsilon } {
		error "relative error $rel_error too large upon evaluating '$get_lst'  ($j / $i)"
	    }
	    if { $rel_error > $rel_max } { set rel_max $rel_error; set maxi $i; set maxj $j }
	    if { $abs_error > $abs_max } { set abs_max $abs_error }
	}
	puts -nonewline "maximum relative deviation upon evaluating '$get_lst': $rel_max  ($maxj / $maxi); maximum absolute deviation: $abs_max "
	if { $absflag==1 } { puts "(zero occured)" } else { puts " " }
    }

    set maxdx 0
    set maxpx 0
    set maxdy 0
    set maxpy 0
    set maxdz 0
    set maxpz 0
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i pr f]
	set tgtF $F($i)
	set dx [expr abs([lindex $resF 0] - [lindex $tgtF 0])]
	set dy [expr abs([lindex $resF 1] - [lindex $tgtF 1])]
	set dz [expr abs([lindex $resF 2] - [lindex $tgtF 2])]

	if { $dx > $maxdx} {
	    set maxdx $dx
	    set maxpx $i
	}
	if { $dy > $maxdy} {
	    set maxdy $dy
	    set maxpy $i
	}
	if { $dz > $maxdz} {
	    set maxdz $dz
	    set maxpz $i
	}
    }
    puts "maximal force deviation in x $maxdx for particle $maxpx, in y $maxdy for particle $maxpy, in z $maxdz for particle $maxpz"
    if { $maxdx > $epsilon || $maxdy > $epsilon || $maxdz > $epsilon } {
	if { $maxdx > $epsilon} {puts "force of particle $maxpx: [part $maxpx pr f] != $F($maxpx)"}
	if { $maxdy > $epsilon} {puts "force of particle $maxpy: [part $maxpy pr f] != $F($maxpy)"}
	if { $maxdz > $epsilon} {puts "force of particle $maxpz: [part $maxpz pr f] != $F($maxpz)"}
	error "force error too large"
    }
} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
