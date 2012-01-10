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
#  Offline Analysis
#

set filename "data/sim_info.dat"
set in [open "$filename" "r"]

while {  [set btitle [ blockfile $in read auto] ] != "eof" } {
 if {$btitle == "variable"} {
     set times [setmd time]
    }
 if { $btitle == "particles"}   {
           # 
           # Here We have an access to all simulation information on given time 
	   # 
     set part0 [part 0 print pos]
     set part0v [part 0 print v]
     puts "time=$times position of p0= $part0 velo = $part0v"
    }
}
close $in
exit
