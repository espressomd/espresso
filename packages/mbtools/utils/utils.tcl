# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
# mbtools::utils -- 
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# Author: Ira
# 

package require ::mmsg 1.0.0
package provide ::mbtools::utils 1.0.0

namespace eval ::mbtools::utils {
    variable verbosity 0
    variable firstconfignum 0

    # The identity number for written out warmup configurations
    variable warmcfg 0

    variable maxnbinteractiontype 0

 


    

}

source [file join [file dirname [info script]] warmup.tcl]
source [file join [file dirname [info script]] topo.tcl]
source [file join [file dirname [info script]] math.tcl]
source [file join [file dirname [info script]] misc.tcl]
source [file join [file dirname [info script]] setup.tcl]



