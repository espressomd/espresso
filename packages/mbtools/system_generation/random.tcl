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
# 
# Routines for creating a random "gas" of lipids
#

namespace eval ::mbtools::system_generation {}

namespace eval ::mbtools::system_generation::random {
    namespace export create_random
}

# ::system::generation::random::create_random --
#
# puts lipids in random positions with random orientation
#
proc ::mbtools::system_generation::random::create_random { args } {
    ::mmsg::send [namespace current] "placing lipids in a random fluid "
    set options {
	{bondl.arg     1.0   "bond length between atoms"  }
	{shuffle "shuffle topology before placement"}
	{exclude.arg "" "a region where no lipids should be placed"}
        {inside.arg "" "a region where lipids should be placed"}
    }
    set usage "Usage: create_random \[bondl:shuffle:exclude:inside]"
    array set params [::cmdline::getoptions args $options $usage]
    
    
    global ::mbtools::system_generation::topology
    global ::mbtools::system_generation::boxl
    
    if { [llength $params(exclude)] != 0 } {
        mmsg::send [namespace current] "Excluding $params(exclude)..."
    }
    
    if { [llength $params(inside)] != 0 } {
        mmsg::send [namespace current] "Inside $params(inside)"
    }
    # First work out how many mol types there are and construct a list
    # of their lengths
    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]

    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    set maxtries 1000

    if { $params(shuffle) } {
        set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }
    
    foreach mol $topology {



	set tries 0
	set tailpos { 0 0 0 }
        set isallowed_ex 0
        set isallowed_in 0
        while { [expr !$isallowed_ex] || [expr !$isallowed_in] } {

            if {  ($tries > $maxtries) } {
                mmsg::err [namespace current] "Could not place molecules: exceeded max number of tries. Please check your conditions."
            }


            # First we choose a random point in space for the tail atom
            lset tailpos 0 [expr $bx*[t_random]]
            lset tailpos 1 [expr $by*[t_random]]
            lset tailpos 2 [expr $bz*[t_random]]
            if { [llength $params(exclude)] != 0 } {            
		# isoutside returns 1 if position is not in the volume
                set isallowed_ex [::mbtools::utils::isoutside $tailpos $params(exclude) ] 
            } else {
                set isallowed_ex 1
            }
            if {[llength $params(inside)] != 0 } {
		# isoutside returns 1 if position is not in the volume
                set isallowed_in ![::mbtools::utils::isoutside $tailpos $params(inside) ] 
            } else {
                set isallowed_in 1
            }
            incr tries

        }

        # Now choose a random orientation vector.  
        lappend orient [expr 2*[t_random]-1]
        lappend orient [expr 2*[t_random]-1]
        lappend orient [expr 2*[t_random]-1] 

        
        ::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
        
        unset tailpos
        unset orient
    }


    # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }

    return

}



