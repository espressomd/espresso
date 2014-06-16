# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#    Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


source "tests_common.tcl"

require_feature "LENNARD_JONES"

puts "----------------------------------------------"
puts "- Testcase analysis.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"

# New configurations for comparison can be added using the WriteNewConfig function
proc WriteNewConfig {in out} {
  global observables listables

  puts -nonewline "Importing checkpoints from '$in' and re-writing them to '$out'... "; flush stdout
  if { 
    [file exists "$in.chk"] } { set chk [open "$in.chk" "r"] 
  } elseif { 
    [file exists "$in"] } { set chk [open "$in" "r"] 
  } else { 
    puts "ERROR: Could not find checkpoint-list $in!\nAborting..."; exit 
  }
  set i 0; 
  set sys_obs [list box_l gamma periodicity skin temperature time time_step]

  while { [eof $chk] == 0 } {
    if { [gets $chk source] > 0 } {
      if { [string compare [lindex [split $source "."] end] "gz"] == 0 } { 
        set f [open "|gzip -cd $source" r]
      } else { 
        set f [open "$source" "r"]
      }
      while { [blockfile $f read auto] != "eof" } {}
      puts -nonewline "."; flush stdout;
      if { [catch { close $f } fid] } {
        puts "Error while closing $file caught: $fid."
      }

      # Get observables from current configuration and store them
      analyze set chains 0 20 30; set observables ""; set listables ""
      
      set mindist1 [analyze mindist];
      set mindist2 [analyze mindist 0 0];
      lappend observables $mindist1 $mindist2
      
      set nbhood [lsort -integer [analyze nbhood 13 2.5]];
      lappend listables $nbhood
      
      set distto [analyze distto 13];
      lappend observables $distto
      
      set energy [analyze energy total];
      lappend observables $energy
      
      set pressure [analyze pressure total];
      lappend observables $pressure
      
      set re [analyze re];
      set re_av [analyze <re>];
      lappend listables $re $re_av
      
      set rg [analyze rg];
      set rg_av [analyze <rg>];
      lappend listables $rg $rg_av

      set rh [analyze rh];
      set rh_av [analyze <rh>];
      lappend listables $rh $rh_av

      set idf [analyze internal_dist];
      set idf_av [analyze <internal_dist>];
      lappend listables $idf $idf_av

      set bdf [analyze bond_dist index 13];
      set bdf_av [analyze <bond_dist> index 13];
      lappend listables $bdf $bdf_av

      set bondl [analyze bond_l];
      set bondl_av [analyze <bond_l>];
      lappend listables $bondl $bondl_av

      set gff [analyze formfactor 1 10 10];
      set gff_av [analyze <formfactor> 1 10 10];
      lappend listables $gff $gff_av

      set g1v [analyze <g1>];
      set g2v [analyze <g2>];
      set g3v [analyze <g3>];
      lappend listables $g1v $g2v $g3v

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
    } 
  }
    close $chk; puts " Done."
}

# Reads the stored configurations from data files
proc ReadConfig { in } {
  global observables listables

  set i 0; 
  set sys_obs [list box_l gamma periodicity skin temperature time time_step]

  set filelist [list "analysis_system.data.00.gz" \
                     "analysis_system.data.01.gz" \
                     "analysis_system.data.02.gz" \
                     "analysis_system.data.03.gz" \
                     "analysis_system.data.04.gz" \
                     "analysis_system.data.05.gz" \
                     "analysis_system.data.06.gz" \
                     "analysis_system.data.07.gz" \
                     "analysis_system.data.08.gz" \
                     "analysis_system.data.09.gz" \
                     "analysis_system.data.10.gz" ]
  
  foreach fname $filelist {
    puts "Reading $fname"
    if { [string compare [ lindex [split $fname "."] end ] "gz" ] ==0 } { 
      set f [open "|gzip -cd $fname" r]
    } else { 
      set f [open "$fname" "r"]
    }
    while { [blockfile $f read auto] != "eof" } {}
    if { [catch { close $f } fid] } {
      puts "Error while closing $f caught: $fid."
    }
  }
}


set epsilon  1e-4
thermostat off
#set tcl_precision  6
set slow     0


# Uncomment to create a new configuration:
  # WriteNewConfig analysis_system.data analysis_system.data2
  # exit 0


###                     ###
#   The analysis script   #
###                     ###

test_catch {
  ReadConfig analysis_system.data; puts "Done."

  analyze set chains 0 20 30\
  set volume [expr pow([lindex [setmd box_l] 0],3)]
  
  set mindist1 "analyze mindist"
  set mindist2 "analyze mindist 0 0"
  lappend get_observables $mindist1 $mindist2

  set nbhood "lsort -integer \[analyze nbhood 13 2.5\]"
  lappend get_listables $nbhood
  
  set distto "analyze distto 13"
  lappend get_observables $distto

  set energy "analyze energy total"
  lappend get_observables $energy
  
  set pressure "analyze pressure total"
  lappend get_observables $pressure
  
  set re "analyze re"
  set re_av "analyze <re>"
  lappend get_listables $re $re_av
  
  set rg "analyze rg"
  set rg_av "analyze <rg>"
  lappend get_listables $rg $rg_av
  
  set rh "analyze rh"
  set rh_av "analyze <rh>"
  lappend get_listables $rh $rh_av

  set idf "analyze internal_dist"
  set idf_av "analyze <internal_dist>"
  lappend get_listables $idf $idf_av

  set bdf "analyze bond_dist index 13"
  set bdf_av "analyze <bond_dist> index 13"
  lappend get_listables $bdf $bdf_av

  set bondl "analyze bond_l"
  set bondl_av "analyze <bond_l>"
  lappend get_listables $bondl $bondl_av

  set gff "analyze formfactor 1 10 10"
  set gff_av "analyze <formfactor> 1 10 10"
  lappend get_listables $gff $gff_av

  set g1v "analyze <g1>"
  set g2v "analyze <g2>"
  set g3v "analyze <g3>"
  lappend get_listables $g1v $g2v $g3v

  for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    set F($i) [part $i pr f]
  }

  integrate 0 recalc_forces

  foreach obs $observables get_obs $get_observables {
    set rel_error [expr abs(([eval $get_obs] - $obs)/$obs)]
    puts "relative deviations upon evaluating '$get_obs': $rel_error  ([eval $get_obs] / $obs)"
    if { $rel_error > $epsilon } {
      error "relative error $rel_error too large upon evaluating '$get_obs'  ([eval $get_obs] / $obs)"
    }
  }

  # Checking stress_tensor and local_stress_tensor calculations
  # the pressure and the stress_tensor may give different answers when the ROTATION option is activated
  # 'analyze stress_tensor' should effectively running on any node

  for { set i 0 } { $i < [setmd n_part] } { incr i } { lappend plist $i }

  #dived box into 4 cuboids - analyse local_stress_tensor in each cuboid  (many bins in each cuboid)
  set box_l [setmd box_l]
  set startx 0.93
  set widthx 0.49
  set starty 0.33
  set widthy 0.57
  set startz 0.10
  set widthz 0.47

  set range_start1 "[expr [lindex $box_l 0]*$startx] 99 99"
  set periodic1 "0 1 1"
  set bins1 "12 1 13"
  set volume1 [expr [lindex $box_l 0]*$widthx * [lindex $box_l 1] * [lindex $box_l 2]]
  set range1 "[expr [lindex $box_l 0]*$widthx] 0 0"
  
  set range_start2 "[expr [lindex $box_l 0]*($startx+$widthx)] [expr [lindex $box_l 1]*$starty] 99"
  set periodic2 "0 0 1"
  set bins2 "1 1 1"
  set range2 "[expr [lindex $box_l 0]*(1-$widthx)] [expr [lindex $box_l 1]* $widthy] 99"
  set volume2 [expr [lindex $box_l 0]*(1-$widthx) * [lindex $box_l 1]*$widthy * [lindex $box_l 2]]
  
  set range_start3 "[expr [lindex $box_l 0]*($startx+$widthx)] [expr [lindex $box_l 1]*($starty+$widthy)] [expr [lindex $box_l 2]*$startz]"
  set periodic3 "0 0 0"
  set bins3 "2 4 20"
  set range3 "[expr [lindex $box_l 0]*(1-$widthx)] [expr [lindex $box_l 1]* (1-$widthy)] [expr [lindex $box_l 2]* $widthz]"
  set volume3 [expr [lindex $box_l 0]*(1-$widthx) * [lindex $box_l 1]*(1-$widthy) * [lindex $box_l 2]*$widthz]
  
  set range_start4 "[expr [lindex $box_l 0]*($startx+$widthx)] [expr [lindex $box_l 1]*($starty+$widthy)] [expr [lindex $box_l 2]*($startz+$widthz)]"
  set periodic4 "0 0 0"
  set range4 "[expr [lindex $box_l 0]*(1-$widthx)] [expr [lindex $box_l 1]* (1-$widthy)] [expr [lindex $box_l 2]* (1-$widthz)]"
  set bins4 " 13 11 3"
  set volume4 [expr [lindex $box_l 0]*(1-$widthx) * [lindex $box_l 1]*(1-$widthy) * [lindex $box_l 2]*(1-$widthz)]

  set local_p_tensor1 [analyze local_stress_tensor \
                        [lindex $periodic1 0] [lindex $periodic1 1] [lindex $periodic1 2]          \
                        [lindex $range_start1 0] [lindex $range_start1 1] [lindex $range_start1 2] \
                        [lindex $range1 0] [lindex $range1 1] [lindex $range1 2]                   \
                        [lindex $bins1 0] [lindex $bins1 1] [lindex $bins1 2] ]
  set local_p_tensor2 [analyze local_stress_tensor \
                        [lindex $periodic2 0] [lindex $periodic2 1] [lindex $periodic2 2]          \
                        [lindex $range_start2 0] [lindex $range_start2 1] [lindex $range_start2 2] \
                        [lindex $range2 0] [lindex $range2 1] [lindex $range2 2]                   \
                        [lindex $bins2 0] [lindex $bins2 1] [lindex $bins2 2] ]
  set local_p_tensor3 [analyze local_stress_tensor \
                        [lindex $periodic3 0] [lindex $periodic3 1] [lindex $periodic3 2]          \
                        [lindex $range_start3 0] [lindex $range_start3 1] [lindex $range_start3 2] \
                        [lindex $range3 0] [lindex $range3 1] [lindex $range3 2]                   \
                        [lindex $bins3 0] [lindex $bins3 1] [lindex $bins3 2] ]     
  set local_p_tensor4 [analyze local_stress_tensor \
                        [lindex $periodic4 0] [lindex $periodic4 1] [lindex $periodic4 2]          \
                        [lindex $range_start4 0] [lindex $range_start4 1] [lindex $range_start4 2] \
                        [lindex $range4 0] [lindex $range4 1] [lindex $range4 2]                   \
                        [lindex $bins4 0] [lindex $bins4 1] [lindex $bins4 2] ]

  set local_p_tensor_sum1 0
  set local_p_tensor_sum2 0
  set local_p_tensor_sum3 0
  set local_p_tensor_sum4 0

  for {set i 1} {$i < [expr [lindex $bins1 0]*[lindex $bins1 1]*[lindex $bins1 2]+1]} {incr i} {
    set local_p_tensor_sum1 [expr $local_p_tensor_sum1               \
                                  +([lindex $local_p_tensor1 $i 1 0] \
                                  + [lindex $local_p_tensor1 $i 1 4] \
                                  + [lindex $local_p_tensor1 $i 1 8])\
                                  *$volume1/[lindex $bins1 0] / [lindex $bins1 1] / [lindex $bins1 2] ]
  }
  for {set i 1} {$i < [expr [lindex $bins2 0]*[lindex $bins2 1]*[lindex $bins2 2]+1]} {incr i} {
    set local_p_tensor_sum2 [expr $local_p_tensor_sum2               \
                                  +([lindex $local_p_tensor2 $i 1 0] \
                                  + [lindex $local_p_tensor2 $i 1 4] \
                                  + [lindex $local_p_tensor2 $i 1 8])\
                                  *$volume2/[lindex $bins2 0] / [lindex $bins2 1] / [lindex $bins2 2] ]
  }
  for {set i 1} {$i < [expr [lindex $bins3 0]*[lindex $bins3 1]*[lindex $bins3 2]+1]} {incr i} {
    set local_p_tensor_sum3 [expr $local_p_tensor_sum3               \
                                  +([lindex $local_p_tensor3 $i 1 0] \
                                  + [lindex $local_p_tensor3 $i 1 4] \
                                  + [lindex $local_p_tensor3 $i 1 8])\
                                  *$volume3/[lindex $bins3 0] / [lindex $bins3 1] / [lindex $bins3 2] ]
  }
  for {set i 1} {$i < [expr [lindex $bins4 0]*[lindex $bins4 1]*[lindex $bins4 2]+1]} {incr i} {
    set local_p_tensor_sum4 [expr $local_p_tensor_sum4               \
                                  +([lindex $local_p_tensor4 $i 1 0] \
                                  + [lindex $local_p_tensor4 $i 1 4] \
                                  + [lindex $local_p_tensor4 $i 1 8])\
                                  *$volume4/[lindex $bins4 0] / [lindex $bins4 1] / [lindex $bins4 2] ]
  }

  set local_p_tensor_sum [expr   $local_p_tensor_sum1 \
                               + $local_p_tensor_sum2 \
                               + $local_p_tensor_sum3 \
                               + $local_p_tensor_sum4]

  set local_p_tensor_total [expr $local_p_tensor_sum/3/[lindex $box_l 0]/[lindex $box_l 1]/[lindex $box_l 2]]

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
    set p_tensorT [expr $p_tensorT + [lindex $total [expr 1 + 4*$i]]/3.0]
    set p_tensorI [expr $p_tensorI + [lindex $ideal [expr 1 + 4*$i]]/3.0]
    set p_tensorF [expr $p_tensorF + [lindex $fene  [expr 2 + 4*$i]]/3.0]
    set p_tensorJ [expr $p_tensorJ + [lindex $lj    [expr 3 + 4*$i]]/3.0]
    set p_tensorJTINTRA [expr $p_tensorJTINTRA + [lindex $nb_tintra    [expr 1 + 4*$i]]/3.0]
    set p_tensorJTINTER [expr $p_tensorJTINTER + [lindex $nb_tinter    [expr 1 + 4*$i]]/3.0]
    set p_tensorJINTRA [expr $p_tensorJINTRA + [lindex $nb_intra    [expr 3 + 4*$i]]/3.0]
  }

  set p_tensor1 $p_tensorT; 
  set rel_error [expr abs(($p_tensor1 - $local_p_tensor_total)/$p_tensor1)]

  puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
        to 'analyze local_stress_tensor...': $rel_error  ($p_tensor1 / $local_p_tensor_total)"
  if { $rel_error > $epsilon } { 
    error "relative error $rel_error too large upon comparing the stress tensors"
  }

  if { ! [regexp "ROTATION" [code_info]]} {
    set p_tot [analyze pressure total]
    set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
    puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
          to 'analyze pressure total': $rel_error  ($p_tot / $p_tensor1)"
    
    if { $rel_error > $epsilon } { 
      error "relative error $rel_error too large upon comparing the pressures" 
    }
    puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
          to 'analyze pressure total': $rel_error  ($p_tot / $p_tensor1)"
    if { $rel_error > $epsilon } {
      error "relative error $rel_error too large upon comparing the pressures" 
    }

    set p_tot [analyze pressure ideal]
    set p_tensor1 $p_tensorI
    set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]

    puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
          to 'analyze pressure ideal': $rel_error  ($p_tot / $p_tensor1)"
    if { $rel_error > $epsilon } { 
      error "relative error $rel_error too large upon comparing the pressures" }
  } else {
    puts "ROTATION is compiled in so total and ideal components\
          of pressure cannot be compared with stress_tensor"
  }
  
  set p_tot [analyze pressure bonded 0]
  set p_tensor1 $p_tensorF
  set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
  puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
        to 'analyze pressure bonded 0': $rel_error  ($p_tot / $p_tensor1)"
  if { $rel_error > $epsilon } { 
    error "relative error $rel_error too large upon comparing the pressures" 
  }
  
  set p_tot [analyze pressure nonbonded 0 0]
  set p_tensor1 $p_tensorJ
  set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
  puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
        to 'analyze pressure nonbonded 0 0': $rel_error  ($p_tot / $p_tensor1)"
  if { $rel_error > $epsilon } { 
    error "relative error $rel_error too large upon comparing the pressures"
  }

  set p_tot [analyze pressure tot_nb_intra]
  set p_tensor1 $p_tensorJTINTRA
  set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
  puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
        to 'analyze pressure tot_nonbonded_intra': $rel_error  ($p_tot / $p_tensor1)"
  if { $rel_error > $epsilon } { 
    error "relative error $rel_error too large upon comparing the pressures"
  }
  
  set p_tot [analyze pressure tot_nb_inter]
  if { $p_tot > 0} {
    set p_tensor1 $p_tensorJTINTER
    set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
    puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
          to 'analyze pressure tot_nonbonded_inter': $rel_error  ($p_tot / $p_tensor1)"
    if { $rel_error > $epsilon } { 
      error "relative error $rel_error too large upon comparing the pressures"
    }
  }

  set p_tot [analyze pressure nb_intra 0 0]
  set p_tensor1 $p_tensorJINTRA
  set rel_error [expr abs(($p_tensor1 - $p_tot)/$p_tot)]
  puts "relative deviations upon comparing trace of 'analyze stress_tensor'\
        to 'analyze pressure nonbonded_intra 0 0': $rel_error  ($p_tot / $p_tensor1)"
  if { $rel_error > $epsilon } { 
    error "relative error $rel_error too large upon comparing the pressures" 
  }
    
  foreach lst $listables get_lst $get_listables {
	  set rel_max 0; set abs_max 0; set absflag 0; set maxi "-"; set maxj "-"
    foreach i $lst j [eval $get_lst] {
      if { [string first "analyze formfactor" "$get_lst"] == 0 || [string first "analyze <formfactor>" "$get_lst"] == 0 } { 
        if { [expr [lindex $i 0]-[lindex $j 0]] < $epsilon } { 
          set i [lindex $i 1]; set j [lindex $j 1] 
        } else { 
          error "different x-coordinates upon comparing '$get_lst'" } 
        }
        
        if { $i!=0 && $j!=0 } { set rel_error [expr abs(($j - $i)/$i)] 
        } else { set rel_error -1; set absflag 1 }

        set abs_error [expr abs($i-$j)]
        
        if { $rel_error > $epsilon } {
          error "relative error $rel_error too large upon evaluating '$get_lst'  ($j / $i)"
        }
        if { $rel_error > $rel_max } { 
          set rel_max $rel_error; set maxi $i; set maxj $j 
        }
        if { $abs_error > $abs_max } { 
          set abs_max $abs_error 
        }
	  }
	  puts -nonewline "maximum relative deviation upon evaluating\
                     '$get_lst': $rel_max  ($maxj / $maxi); maximum\
                     absolute deviation: $abs_max "

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
  puts "maximal force deviation in x $maxdx for particle $maxpx,\
        in y $maxdy for particle $maxpy, in z $maxdz for particle $maxpz"
  if { $maxdx > $epsilon || $maxdy > $epsilon || $maxdz > $epsilon } {
	if { $maxdx > $epsilon} {puts "force of particle $maxpx: [part $maxpx pr f] != $F($maxpx)"}
	if { $maxdy > $epsilon} {puts "force of particle $maxpy: [part $maxpy pr f] != $F($maxpy)"}
	if { $maxdz > $epsilon} {puts "force of particle $maxpz: [part $maxpz pr f] != $F($maxpz)"}
	error "force error too large"
    }
}

exit 0
