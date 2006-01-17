#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
# 
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 
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
	set nbhood [lsort -integer [analyze nbhood 13 2.5]]; lappend listables $nbhood
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
	set f [open "|gzip - > $out.[eval format %02d $i].gz" w]
	blockfile $f write variable $sys_obs
	blockfile $f write tclvariable observables
	blockfile $f write tclvariable listables
	blockfile $f write interactions
	blockfile $f write integrate
	blockfile $f write thermostat
	blockfile $f write particles "id pos type v f"
	blockfile $f write bonds
	blockfile $f write configs 1
	close $f
	incr i
    } }
    close $chk; puts " Done."
}

set epsilon  1e-4
setmd temp   0.0
setmd gamma  0.0
#set tcl_precision  6
set slow     0

# rewrite analysis_system.data analysis_system.data2; exit 0

if { [catch {
    puts -nonewline "Reading the checkpoint... "; flush stdout
    checkpoint_read analysis_system.data; puts " Done."

    analyze set chains 0 20 30; set volume [expr pow([lindex [setmd box_l] 0],3)]
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
	if { [string first "analyze pressure" $get_obs]==0 } { if { [regexp "ROTATION" [code_info]]} { 
	    puts -nonewline "ROTATION compiled in => adjusting stored pressure $obs (f=3) by current ideal one ([analyze pressure ideal]) "
	    set obs [expr $obs - [analyze pressure ideal]]; puts "to $obs (f=6)"
	} }
	set rel_error [expr abs(([eval $get_obs] - $obs)/$obs)]
	puts "relative deviations upon evaluating '$get_obs': $rel_error  ([eval $get_obs] / $obs)"
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error too large upon evaluating '$get_obs'  ([eval $get_obs] / $obs)"
	}
    }
    
    # Checking p_IK1 stress tensor calculation
    if { [setmd n_nodes]==1 || $slow==1 } {
	# since 'analyze p_IK1' is effectively running on the master node only, it will be really slow for >1 nodes;
	# hence it is not really necessary to check it there - but if you like, set 'slow' to 1
	for { set i 0 } { $i < [setmd n_part] } { incr i } { lappend plist $i }
	set p_IK1_full [analyze p_IK1 $volume $plist 0];
	set total [lindex $p_IK1_full 0]
	set ideal [lindex $p_IK1_full 1]
	set fene  [lindex $p_IK1_full 2]
	set lj    [lindex $p_IK1_full 3]
        set p_IKT 0; set p_IKI 0; set p_IKF 0; set p_IKJ 0
	for {set i 0} {$i < 3} {incr i} {
	    set p_IKT [expr $p_IKT + [lindex $total [expr 1 + 4*$i]]]
	    set p_IKI [expr $p_IKI + [lindex $ideal [expr 1 + 4*$i]]]
	    set p_IKF [expr $p_IKF + [lindex $fene  [expr 2 + 4*$i]]]
	    set p_IKJ [expr $p_IKJ + [lindex $lj    [expr 3 + 4*$i]]]
	}
        set p_tot [analyze pressure total]; set p_IK1 $p_IKT; set rel_error [expr abs(($p_IK1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze p_IK1' to 'analyze pressure total': $rel_error  ($p_tot / $p_IK1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure ideal]; set p_IK1 $p_IKI; set rel_error [expr abs(($p_IK1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze p_IK1' to 'analyze pressure ideal': $rel_error  ($p_tot / $p_IK1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure fene 0]; set p_IK1 $p_IKF; set rel_error [expr abs(($p_IK1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze p_IK1' to 'analyze pressure fene 0': $rel_error  ($p_tot / $p_IK1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure lj 0 0]; set p_IK1 $p_IKJ; set rel_error [expr abs(($p_IK1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze p_IK1' to 'analyze pressure lj 0 0': $rel_error  ($p_tot / $p_IK1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
    }

    # Checking stress_tensor calculation
    # 'analyze stress_tensor' should effectively running on any node
    #if { [setmd n_nodes]==1 || $slow==1 } {
	for { set i 0 } { $i < [setmd n_part] } { incr i } { lappend plist $i }
	set p_tensor_full [analyze stress_tensor];
	set total [lindex $p_tensor_full 0]
	set ideal [lindex $p_tensor_full 1]
	set fene  [lindex $p_tensor_full 2]
	set lj    [lindex $p_tensor_full 3]
	set nb_tintra    [lindex $p_tensor_full 4]
	set nb_tinter    [lindex $p_tensor_full 5]
	set nb_intra    [lindex $p_tensor_full 6]
        set p_tensorT 0; set p_tensorI 0; set p_tensorF 0; set p_tensorJ 0; 
	# checking intra- and inter- molecular non-bonded contribution to stress tensor
	set p_tensorJTINTRA 0; set p_tensorJTINTER 0; set p_tensorJINTRA 0; 
	for {set i 0} {$i < 3} {incr i} {
	    set p_tensorT [expr $p_tensorT + [lindex $total [expr 1 + 4*$i]]]
	    set p_tensorI [expr $p_tensorI + [lindex $ideal [expr 1 + 4*$i]]]
	    set p_tensorF [expr $p_tensorF + [lindex $fene  [expr 2 + 4*$i]]]
	    set p_tensorJ [expr $p_tensorJ + [lindex $lj    [expr 3 + 4*$i]]]
	    set p_tensorJTINTRA [expr $p_tensorJTINTRA + [lindex $nb_tintra    [expr 1 + 4*$i]]]
	    set p_tensorJTINTER [expr $p_tensorJTINTER + [lindex $nb_tinter    [expr 1 + 4*$i]]]
	    set p_tensorJINTRA [expr $p_tensorJINTRA + [lindex $nb_intra    [expr 3 + 4*$i]]]
	}
        set p_tot [analyze pressure total]; set p_tensor1 $p_tensorT; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure total': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure ideal]; set p_tensor1 $p_tensorI; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure ideal': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure fene 0]; set p_tensor1 $p_tensorF; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure fene 0': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure lj 0 0]; set p_tensor1 $p_tensorJ; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure lj 0 0': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure tot_nb_intra]; set p_tensor1 $p_tensorJTINTRA; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure tot_nonbonded_intra': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
        set p_tot [analyze pressure tot_nb_inter]; 
	if { $p_tot > 0} {
	  set p_tensor1 $p_tensorJTINTER; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	  puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure tot_nonbonded_inter': $rel_error  ($p_tot / $p_tensor1)"
	  if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
	}
        set p_tot [analyze pressure nb_intra 0 0]; set p_tensor1 $p_tensorJINTRA; set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
	puts "relative deviations upon comparing trace of 'analyze stress_tensor' to 'analyze pressure nonbonded_intra 0 0': $rel_error  ($p_tot / $p_tensor1)"
	if { $rel_error > $epsilon } { error "relative error $rel_error too large upon comparing the pressures" }
    #}
       
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
