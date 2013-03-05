# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#   
# ESPResSo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# This set of routines are used to set constraints using the expresso
# constraint feature
#

namespace eval ::mbtools::system_generation {}

proc ::mbtools::system_generation::create_sphericalconstraint { topo box_l args } { 
    mmsg::send [namespace current] "making a spherical constraint "
    # ---- Process Command Line Args --------------------------- #
    
    set options {
	{center.arg  { 0.0 0.0 0.0 } "location of the sphere relative to the box center" }
	{radius.arg  1.0      "radius of the constraint" }
	{dir.arg    1  "constrain out or in -1 is inward" }
	{type.arg  0  "the type of the constraint for the purposes of setting interactions"}

    }
    set usage "Usage: create_uniform_sphere alipid: "
    array set params [::cmdline::getoptions args $options $usage]

    #Setup the center so it is offset from the box center
    set center [list 0.0 0.0 0.0 ]

    #Construct the real center
    lset center 0 [expr [lindex $params(center) 0] + [lindex $box_l 0]/(2.0)]
    lset center 1 [expr [lindex $params(center) 1] + [lindex $box_l 1]/(2.0)]
    lset center 2 [expr [lindex $params(center) 2] + [lindex $box_l 2]/(2.0)]


    constraint sphere center [lindex $center 0] [lindex $center 1] [lindex $center 2] radius $params(radius) direction $params(dir) type $params(type)


    return
    
    


}


