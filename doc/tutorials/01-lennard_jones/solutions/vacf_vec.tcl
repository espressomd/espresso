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
#  Simple Velocity Auto-Correlation-Function; Direct Method
#  February 2007
#  by m@suzen 
#
set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/vacf.dat" "w"]
set origin 0
puts $out "#"
puts $out "#"
puts $out "# Velocity Autocorrelation Functions"
puts $out "#"
puts  "Reading frame"
set frame_count 0
while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
   if {$btitle == "variable"} {
     set times [setmd time]
     set pmax [setmd max_part]
    }
    if { $btitle == "particles" && $pmax > 0}   {
        puts -nonewline "Reading frame: $frame_count \r"
      # Origin Value
          if {$origin == 0} {
             for {set i 0} {$i < $pmax} {incr i} {
                set v [part $i print v]
                lappend pvx0 [lindex $v 0] ;# initial velocity components
                lappend pvy0 [lindex $v 1]
                lappend pvz0 [lindex $v 2]
              }
	      set origin 1
	    }
        # AFC 
              set acf 0.0
	      set vx [list]
	      set vy [list]
	      set vz [list]
             for {set i 0} {$i < $pmax} {incr i} {
                set v [part $i print v]
                lappend vx [lindex $v 0]
                lappend vy [lindex $v 1]
                lappend vz [lindex $v 2]
              }
              set acf [expr [vecdot_product $pvx0 $vx]+[vecdot_product $pvy0 $vy]+[vecdot_product $pvz0 $vz]]
	      set acf [expr $acf/$pmax]
	      puts $out "$times   $acf  "
             set frame_count [expr $frame_count+1]
       }
}
close $in
close $out
exit
