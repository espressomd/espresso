#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;
# 
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
# 
set nrep { 1000 1000 1000 1000 1000 1000 1000 1000 }

proc blubb {vec} {
   return [expr -1 * log([lindex $vec 1] / [lindex $vec 0])]
}

puts "-------------------------------------------"
puts "- Testcase uwerr.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------"
set df [open "uwerr_test.data" r]
while {![eof $df]} {
   gets $df row
   lappend data [split $row " "]
}
close $df
set data [lrange $data 0 end-1]

puts "Expected values:"
puts "0.190161129416 0.0149872743495 0.00120248945994 8.70337780606 1.27314416767 0.992579964046"
puts [uwerr $data $nrep blubb plot]
