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
#  Offline Analysis for Radial Distribution Function
#
#############################################################
#                                                           #
#  Determines the radial distribution function              #
#  from the data in sim_info.dat file                       #
#                                                           # 
#############################################################

# set the input and output streams
set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/rdf.dat" "w"]

# here we set up list in which we are going to bin the 
# distance dependent rdf results, we set the number of 
# bins to 100 here and set the value zero in each bin
set rbin 100
for {set i 0} {$i < $rbin} {incr i} { 
  lappend avg_rdf 0; 
}

# set the type of particles we want the RDF for and
# the lower boundary for the radial regime
set group 0
set rmin 0.0

# read in all the data in file and set the counter
set frame_no 1; 
while { [set btitle [ blockfile $in read auto ] ] != "eof" } {

  if {$btitle == "variable"} {

    # here we determine the simulation time and the box length,
    # note that we implicitly assume that the box is cubic. The
    # maximum radius to which we can determine the RDF is half
    # the box size. Otherwise we run into periodicity problems
    set times [setmd time]
    set boxl [lindex [setmd box_l] 0]
    set rmax [expr $boxl/2.0]

    puts -nonewline "frame $frame_no\r"
    flush stdout
  }

  if { $btitle == "particles"} {

    # in the particles block we can invoke the analyze rdf 
    # command. We want to analyze the distribution around 
    # particles of type $group of particles with type $group
    # (e.g., in charged systems one can obtain multiple RDFs)
    # we want to do this in the radial interval [$rmin,$rmax]
    # with an number of equidistant bins equal to $rbin.
    set drdf [analyze rdf $group $group $rmin $rmax $rbin ]

    # we now extract the RDF from the output of the analyze 
    # routine. We also determine the size of the array.
    set data_rdf [lindex $drdf 1]
    set dsize [llength $data_rdf]

    # we create two empty lists, the first will contain the
    # radial coordinates, the second the RDF values
    set rlist [list]
    set rdflist [list]

    # here we fill the lists we created earlier. Note that
    # the algorithm could be more efficient, since the first
    # would only need to be constructed once
    for {set i 0} {$i <$dsize} {incr i} {
      lappend  rlist [lindex $data_rdf $i 0]
      lappend  rdflist [lindex $data_rdf $i 1]
    }

    # to determine the average, we add the values in the
    # rdflist to a running total, making use of the 
    # vecadd function, which allows for element-wise addition
    # of arrays of the same length
    set avg_rdf [vecadd $avg_rdf $rdflist]

    # finally we increment the frame number, again this must
    # be done within either of the if statements, otherwise
    # the number of frames is overcounted by a factor of 9 ! 
    incr frame_no ;
  }
}

# needed to show the last frame number
puts "\n"

# here the element-wise vecscale is used to divide the 
# RDF totals by the total number of frames
set avg_rdf [vecscale [expr 1.0/$frame_no]  $avg_rdf]

# finally the information is printed to the file. Note
# that here the foreach command may be used since the 
# two arrays are of the same length
foreach r $rlist value $avg_rdf { puts $out "$r $value" }

# close the streams and exit ESPResSo
close $in
close $out
exit
