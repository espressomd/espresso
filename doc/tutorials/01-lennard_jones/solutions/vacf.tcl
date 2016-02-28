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
#  Simple Velocity Auto-Correlation-Function; Direct Method
#
set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/vacf.dat" "w"]
set origin 0
puts $out "#"
puts $out "#"
puts $out "# Velocity Autocorrelation Functions"
puts $out "#"
set frame_count 0
while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
   if {$btitle == "variable"} {
     set times [setmd time]
     set pmax [setmd max_part]
    }
    if { $btitle == "particles" && $pmax > 0}   {
        puts "Reading frame: $frame_count "
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
             for {set i 0} {$i < $pmax} {incr i} {
                set v [part $i print v]
                set vx [lindex $v 0]
                set vy [lindex $v 1]
                set vz [lindex $v 2]
		set vx0 [lindex $pvx0 $i]
		set vy0 [lindex $pvy0 $i]
		set vz0 [lindex $pvz0 $i]
                set acf [expr $vx*$vx0+$vy*$vy0+$vz*$vz0+$acf]
              }
	      set acf [expr $acf/$pmax]
	      puts $out "$times   $acf  "
             set frame_count [expr $frame_count+1]
       }
}
close $in
close $out
exit
