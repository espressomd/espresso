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

set filename "data/sim_info.dat"
set in [open "$filename" "r"]
set out [open "data/rdf.dat" "w"]

set rbin 100
for {set i 0} {$i < $rbin} {incr i} { 
  lappend avg_rdf 0; 
}

set group 0
set rmin 0.0

set frame_no 1; 
while { [set btitle [ blockfile $in read auto ] ] != "eof" } {

  if {$btitle == "variable"} {

    set times [setmd time]
    set boxl [lindex [setmd box_l] 0]
    set rmax [expr $boxl/2.0]

    puts -nonewline "frame $frame_no\r"
    flush stdout
  }

  if { $btitle == "particles"} {

    set drdf [analyze rdf $group $group $rmin $rmax $rbin ]

    set data_rdf [lindex $drdf 1]
    set dsize [llength $data_rdf]

    set rlist [list]
    set rdflist [list]

    for {set i 0} {$i <$dsize} {incr i} {
      lappend  rlist [lindex $data_rdf $i 0]
      lappend  rdflist [lindex $data_rdf $i 1]
    }

    set avg_rdf [vecadd $avg_rdf $rdflist]

    incr frame_no ;
  }
}

puts "\n"

set avg_rdf [vecscale [expr 1.0/$frame_no]  $avg_rdf]

foreach r $rlist value $avg_rdf { puts $out "$r $value" }

close $in
close $out
exit
