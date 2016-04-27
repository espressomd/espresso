# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#
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
