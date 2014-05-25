# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "FFTW"

puts "---------------------------------------------------------------"
puts "- Testcase p3m_simple_noncubic.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 10 10 10

setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m 3 32 5 0.001
# in k-space:
inter coulomb 1.0 p3m 0.0 32 5 1.0

integrate 0
puts [ part 0 print f ]
set f1 [ lindex [part 0 print f] 0 ]

setmd box_l 20 10 10

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
part 2 pos 14 5 5 q +1 v 0 0 0
part 3 pos 16 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m 3 32 5 0.001
# in k-space:

if { [catch {

    inter coulomb 1.0 p3m 0.0 64 5 1.0
    integrate 0

    puts [ part 0 print f ]
    set f2 [ lindex [part 0 print f] 0 ]
    
    if { abs( $f1 - $f2 )/$f1 < 1e-5 } {
	puts OK
    } else {
	error "P3M noncubic test failed"
    }

} res ] } {
    error_exit $res
}

exit 0
