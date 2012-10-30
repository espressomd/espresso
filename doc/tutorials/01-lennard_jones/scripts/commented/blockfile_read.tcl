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
#  Offline Analysis
#
#############################################################
#                                                           #
#  Reads in data saved by the save_sim in lj_tutorial.tcl   #
#                                                           # 
#############################################################

# set the filename and open a read stream to it
set filename "data/sim_info.dat"
set in [open "$filename" "r"]

# read in the data blocks from the block file
while { [set btitle [ blockfile $in read auto ] ] != "eof" } {

  # set the time. The if statement prevents the time
  # to be set in the other blocks. That is, it is only 
  # set once per snapshot, instead of 9 times
  if {$btitle == "variable"} {
    set times [setmd time]
  }

  # set velocity of the first particle, 
  # which is in "particles" block, so it can only
  # be set after that is read in
  if { $btitle == "particles"} {
    set part0 [part 0 print pos]
    set part0v [part 0 print v]
    puts "time = $times"
    puts "position of the first particle = $part0"
    puts "velocity of the first particle = $part0v\n"
  }
}

# close the stream and exit ESPResSo
close $in
exit
