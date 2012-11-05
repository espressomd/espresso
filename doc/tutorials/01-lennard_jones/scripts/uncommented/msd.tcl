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
#  Offline Analysis for Mean Square Deviation
#
#############################################################
#                                                           #
#  Determines the mean square displacement                  #
#  from the data in sim_info.dat file                       #
#                                                           # 
#############################################################

set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/msd.dat" "w"]

set origin 0

set frame_no 0
while { [set btitle [ blockfile $in read auto] ] != "eof" } {

  set msd 0.0
  if { $btitle == "particles"} {

    if {$origin == 0} {
      set time0 [setmd time]
    }
    set times [setmd time]

    set pmax [setmd n_part]

    if {$origin == 0} {
      for {set i 0} {$i < $pmax} {incr i} {
        set pos [part $i print pos]
        lappend vec_r0 $pos
      }

      set origin 1
    }

    if {$origin == 1} {
      for {set i 0} {$i < $pmax} {incr i} {
        set pos [part $i print pos]
        set pos0 [lindex $vec_r0 $i]
        set delta_r [vecsub $pos $pos0]
        set msd [expr $msd + [veclensqr $delta_r] ]
      }
      set msd [expr $msd/$pmax]
    }

    puts $out "[expr $times - $time0] $msd"

    incr frame_no
    puts -nonewline "frame $frame_no\r"
    flush stdout
  }
}

puts "\n"

close $in
close $out
exit
