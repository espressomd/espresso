#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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
set cnt 0
for {set i 0} {$i < 100} {incr i} { lappend avg_rdf 0 }

foreach filename [lrange $argv 1 end] {

    set f [open $filename "r"]
    while { [blockfile $f read auto] != "eof" } {}
    close $f

    set rdf [analyze rdf 0 1 0.9 [expr $box_l/2] 100]
    set rlist ""
    set rdflist ""
    foreach value [lindex $rdf 1] {
	lappend rlist [lindex $value 0]
	lappend rdflist [lindex $value 1]
    }
    
    set avg_rdf [vecadd $avg_rdf $rdflist]
    incr cnt
}

set avg_rdf [vecscale [expr 1.0/$cnt] $avg_rdf]

set plot [open "rdf.data" "w"]
puts $plot "\# r rdf(r)"
foreach r $rlist rdf $avg_rdf { puts $plot "$r $rdf" }
close $plot
