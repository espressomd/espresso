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

# set the input and output streams
set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/msd.dat" "w"]

# required to extract the first point in the simulation
# i.e. form which we calculate all time differences
set origin 0

# read in all the data in file and set the counter
set frame_no 0
while { [set btitle [ blockfile $in read auto] ] != "eof" } {

  set msd 0.0
  if { $btitle == "particles"} {

    # set the time and catch the first time to fix 
    # the origin of the time axis; the production run
    # does not start at t=0 !
    if {$origin == 0} {
      set time0 [setmd time]
    }
    set times [setmd time]

    # The number of particles can only be set when
    # the particles have been read in
    set pmax [setmd n_part]

    # set the values for the posititions
    # for the 'origin' point of interest
    if {$origin == 0} {
      for {set i 0} {$i < $pmax} {incr i} {
        set pos [part $i print pos]
        lappend vec_r0 $pos
      }

      # we have the origin point, now let us 
      # consider the other points
      set origin 1
    }

    # now fetch the current value of the 
    # position vectors, subtract, square and
    # average this value over all particles
    if {$origin == 1} {
      for {set i 0} {$i < $pmax} {incr i} {
        set pos [part $i print pos]
        set pos0 [lindex $vec_r0 $i]
        set delta_r [vecsub $pos $pos0]
        set msd [expr $msd + [veclensqr $delta_r] ]
      }
      set msd [expr $msd/$pmax]
    }

    # output the msd point we are currently working on
    puts $out "[expr $times - $time0] $msd"

    # show the progress
    incr frame_no
    puts -nonewline "frame $frame_no\r"
    flush stdout
  }
}

# needed to show the last frame number
puts "\n"

# close the streams and exit ESPResSo
close $in
close $out
exit
