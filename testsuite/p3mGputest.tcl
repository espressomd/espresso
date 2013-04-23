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



puts "---------------------------------------------------------------"
puts "- Testcase p3m_simple_noncubic.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 10 10 10

#lbfluid gpu dens 1.0 visc 1.0 agrid 1 friction 0.00000000001 tau 0.01
setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m gpu 3 32 5 0.001
# in k-space:
#puts [inter coulomb 1.0 p3m gpu tunev2 accuracy 0.001]
inter coulomb 1.0 p3m gpu 3.0 32 5 1.0

integrate 0
puts "after integ 0 forces are:"
puts [ part 0 print f ]
puts [ part 1 print f ]
integrate 1
puts "after integ 1 forces are:"
puts [ part 0 print f ]
puts [ part 1 print f ]

set f1 [ lindex [part 0 print f] 0 ]

invalidate_system 

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0

exit 0
