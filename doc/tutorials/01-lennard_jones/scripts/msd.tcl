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
#  Offline Analysis for Mean Square Deviation
#

set filename "data/sim_info.dat"
set in [open "$filename" "r"]

set out [open "data/msd.dat" "w"]
set origin 0
puts $out "#"
puts $out "#"
puts $out "# Mean Square Deviation"
puts $out "#"

set frame_no -1;
while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
 if {$btitle == "variable"} {
     set times [setmd time]
     set pmax [setmd max_part]
     incr frame_no
    }
 set msd 0.0;
 if { $btitle == "particles"}   {
  puts "frame_no $frame_no pmax $pmax"
         if {$frame_no == 1} {
	     for {set i 0} {$i < $pmax} {incr i} {
	       set pos [part $i print pos]
               lappend vec_r0 $pos
	      }
	  }
         if {$frame_no > 1} {
	     for {set i 0} {$i < $pmax} {incr i} {
	       set pos [part $i print pos]
	       set pos0 [lindex $vec_r0 $i]
	       set delta_r [vecsub  $pos $pos0]
	        puts "delta_r $delta_r"
	       set msd [expr $msd + [veclensqr $delta_r] ]
	      }
	       set msd [expr $msd/$pmax]
	       puts $out "$frame_no  $msd"
	  }
    }
}
close $in
close $out
exit
