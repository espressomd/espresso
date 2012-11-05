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

# set the input and output streams
set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/vacf.dat" "w"]

# required to extract the first point in the simulation
# i.e. form which we calculate all time differences
set origin 0

# read in all the data in file and set the counter
set frame_count 0
while {  [set btitle [ blockfile $in read auto ] ] != "eof" } {

  # if we have the particles block and there are
  # particles in inside it we do the following
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

    # start with frame 1
    incr frame_count
    # output the frame we are currently working on
    puts -nonewline "frame $frame_count\r"
    flush stdout

    # set the values for the velocities in three
    # lists for the 'origin' point of interest
    if {$origin == 0} {
      for {set i 0} {$i < $pmax} {incr i} {
        set v [part $i print v]
        lappend pvx0 [lindex $v 0]
        lappend pvy0 [lindex $v 1]
        lappend pvz0 [lindex $v 2]
      }

      # we have the origin point, now let us 
      # consider the other points
      set origin 1
    }

    # for each point after the initial point
    # we determine the dot product between 
    # the velocities of a particle at time 
    # t = 0 and time t and average these over
    # all particles
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

    # finally we write this result to disk
    puts $out "[expr $times - $time0] $acf"
  }
}

# needed to show the last frame number
puts "\n"

# close the streams and exit ESPResSo
close $in
close $out
exit
