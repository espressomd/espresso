# Copyright (C) 2010,2011,2012 The ESPResSo project
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
#  Simple Velocity Auto-Correlation-Function; Direct Method
#
#############################################################
#                                                           #
#  Determines the velocity auto-correlation function        #
#  from the data in sim_info.dat file                       #
#                                                           # 
#############################################################

set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/vacf.dat" "w"]

set origin 0

set frame_count 0
while {  [set btitle [ blockfile $in read auto ] ] != "eof" } {

  if { $btitle == "particles"} {

    if {$origin == 0} {
      set time0 [setmd time]
    }
    set times [setmd time]

    set pmax [setmd n_part]

    incr frame_count

    puts -nonewline "frame $frame_count\r"
    flush stdout

    if {$origin == 0} {
      for {set i 0} {$i < $pmax} {incr i} {
        set v [part $i print v]
        lappend pvx0 [lindex $v 0]
        lappend pvy0 [lindex $v 1]
        lappend pvz0 [lindex $v 2]
      }

      set origin 1
    }

    set acf 0.0
    for {set i 0} {$i < $pmax} {incr i} {
      set v [part $i print v]
      set vx [lindex $v 0]
      set vy [lindex $v 1]
      set vz [lindex $v 2]
      set vx0 [lindex $pvx0 $i]
      set vy0 [lindex $pvy0 $i]
      set vz0 [lindex $pvz0 $i]
      set acf [expr $acf + ($vx*$vx0 + $vy*$vy0 + $vz*$vz0)]
    }
    set acf [expr $acf/$pmax]

    puts $out "[expr $times - $time0] $acf"
  }
}

puts "\n"

close $in
close $out
exit
