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
#  Offline Analysis for Radial Distribution Function
#

set filename "data/sim_info.dat"
set in [open "$filename" "r"]

set out [open "data/rdf.dat" "w"]
set origin 0
puts $out "#"
puts $out "#"
puts $out "# Radial Distribution Function "
puts $out "#"

set frame_no 1;
# rdf collect
for {set i 0} {$i < 100} {incr i}  { lappend avg_rdf 0; }
while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
 if {$btitle == "variable"} {
     set times [setmd time]
     set boxl [lindex [setmd box_l] 0]
      puts "frame $frame_no"
    }
 if { $btitle == "particles"}   {
           # 
           # Here We have an access to all simulation information on given time 
	   # 
	      set group 0
	      set rmin 0.0
	      set rmax [expr $boxl/2.0]
	      set rbin 100
	      set drdf [analyze rdf $group $group $rmin $rmax $rbin ]
	      set data_rdf [lindex $drdf 1]
	      set dsize [llength $data_rdf]
	          set rlist [list]
	          set rdflist [list]
              for {set i 0} {$i <$dsize} {incr i} {
	        lappend  rlist [lindex $data_rdf $i 0]
	        lappend  rdflist [lindex $data_rdf $i 1]
	      }
	       set avg_rdf [vecadd $avg_rdf  $rdflist]
       incr frame_no ;
    }
}
set avg_rdf [vecscale [expr 1.0/$frame_no]  $avg_rdf]
 foreach r $rlist value $avg_rdf { puts $out "$r  $value" }
close $in
close $out
exit
