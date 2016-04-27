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
# ::mbtools::analysis::analyze_energy --
#
# Calculate the total energy of the system and break it into components
#   

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::energy {
    variable av_components_en 0
    variable av_total_en 0
    variable av_kin_en 0
    variable av_fene_en 0
    variable av_harm_en 0
    variable av_nb_en 0
    variable av_en_i 0

    variable f_tvsen
    variable verbose

    namespace export printav_energy
    namespace export setup_energy
    namespace export analyze_energy
    namespace export resetav_energy
}

proc ::mbtools::analysis::energy::resetav_energy { } {
    variable av_components_en 
    variable av_total_en 
    variable av_kin_en 
    variable av_fene_en 
    variable av_harm_en 
    variable av_nb_en 
    variable av_en_i 
    for {set i 0 } {$i < [llength $av_components_en] } { incr i } {
	lset av_components_en $i 0.0
    }
    set av_total_en 0
    set av_kin_en 0
    set av_fene_en 0
    set av_harm_en 0
    set av_nb_en 0
    set av_en_i 0
}

proc ::mbtools::analysis::energy::printav_energy { } {
    variable av_components_en 
    variable av_total_en 
    variable av_kin_en 
    variable av_fene_en 
    variable av_harm_en 
    variable av_nb_en 
    variable av_en_i     
    variable f_tvsen
    global ::mbtools::analysis::time
    if { $av_en_i > 0 } {
	puts -nonewline $f_tvsen "$time [expr $av_total_en/(1.0*$av_en_i)] [expr $av_fene_en/(1.0*$av_en_i)] [expr $av_kin_en/(1.0*$av_en_i)] [expr $av_harm_en/(1.0*$av_en_i)] [expr $av_nb_en/(1.0*$av_en_i)]"
	foreach comp $av_components_en {
	    puts -nonewline $f_tvsen " [expr $comp/(1.0*$av_en_i)]"
	}
	puts $f_tvsen ""
    }
    flush $f_tvsen
}

proc ::mbtools::analysis::energy::setup_energy { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsen
    variable verbose
    variable av_components_en

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_energy verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_energy$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    ::mmsg::debug  [namespace current]  "opening $outputdir/time_vs_energy$suffix "
    set f_tvsen [open "$outputdir/time_vs_energy$suffix" $iotype ]
    # First work out the names of components
    set raw [analyze energy]
    
    if { $newfile || $iotype == "w"} {
	puts $f_tvsen "\# Components of the energy "
	puts -nonewline $f_tvsen "\# Time total_en fene_en kinetic_en harmonic_en nonbonded_en"
    }
    unset av_components_en
    for { set k 0 } { $k < [llength $raw ] } { incr k } {
	set tmp [lindex $raw $k]
	set ntmp [llength $tmp]	
	if { [ regexp "nonbonded" $tmp ] } {
	    puts -nonewline $f_tvsen " [lrange $tmp 0 end-1]"
	    lappend av_components_en 0.0
	}
    }
    puts $f_tvsen ""
    flush $f_tvsen
}


proc ::mbtools::analysis::energy::analyze_energy {  } {
    ::mmsg::send [namespace current] "analyzing energy"
    variable av_components_en
    variable av_total_en
    variable av_kin_en
    variable av_fene_en
    variable av_harm_en
    variable av_nb_en
    variable av_en_i
    variable verbose
    set energy_all [analyze energy]
#    puts $energy_all
    set nb_en 0
    set nbcount 0
    for { set i 0 } { $i < [llength $energy_all ] } { incr i } {
	set tmp [lindex $energy_all $i]
	set ntmp [llength $tmp]	
	if { [ regexp "energy" $tmp ] } {
	    set total_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "kinetic" $tmp ] } {
	    set kin_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "FENE" $tmp ] } {
	    set fene_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "HARMONIC" $tmp ] } {
	    set harm_en [lindex $tmp [expr $ntmp -1]]
	}
	if { [ regexp "nonbonded" $tmp ] } {
	    set nb_en [expr $nb_en + [lindex $tmp [expr $ntmp -1]]]
#	    puts "$nbcount [llength $av_components_en]"
	    lset av_components_en $nbcount [expr [lindex $av_components_en $nbcount] + [lindex $tmp [expr $ntmp -1]] ]

	    incr nbcount
	}
	
    }
    incr av_en_i
    set av_total_en  [expr $av_total_en + $total_en/(1.0)]
    set av_kin_en  [expr $av_kin_en + $kin_en/(1.0)]
    set av_fene_en  [expr $av_fene_en + $fene_en/(1.0)]
    set av_harm_en  [expr $av_harm_en + $harm_en/(1.0)]
    set av_nb_en  [expr $av_nb_en + $nb_en/(1.0)]
    if { $verbose } {
	::mmsg::send [namespace current] "energy: [expr $av_total_en/($av_en_i*1.0)] kinetic: [expr $av_kin_en/($av_en_i*1.0)] FENE: [expr $av_fene_en/($av_en_i*1.0)] HARMONIC: [expr $av_harm_en/($av_en_i*1.0)] nonbonded: [expr $av_nb_en/($av_en_i*1.0)]"
    }
    


    ::mmsg::debug [namespace current] "done"
    

}


