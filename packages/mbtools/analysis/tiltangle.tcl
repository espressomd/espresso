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
##############################################################################
#                                                                            #   
#Measure the angle between the lipids around a protein and the bilayer normal#
#                                                                            #  
#Author: Gregoria                                                            #
#                                                                            #
##############################################################################

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::tiltangle {
    variable av_toptilt { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 }
    variable av_tilt_i 0
    variable f_tvstilt
    variable verbose

    namespace export setup_tiltangle
    namespace export analyze_tiltangle
    namespace export printav_tiltangle
    namespace export resetav_tiltangle
}

proc ::mbtools::analysis::resetav_tiltangle { } {
    variable toptilt
    variable av_tilt_i
    set av_toptilt {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0}
    set av_tilt_i 0
}

proc ::mbtools::analysis::printav_tiltangle { } { 
    variable av_toptilt
    variable av_tilt_i
    variable f_tvstilt
    global ::mbtools::analysis::time

    if { $av_tilt_i > 0 } {
                   
	set avtoptilt0 [expr [lindex $av_toptilt 0]/($av_tilt_i*1.0)] 
	set avtoptilt1 [expr [lindex $av_toptilt 1]/($av_tilt_i*1.0)]
	set avtoptilt2 [expr [lindex $av_toptilt 2]/($av_tilt_i*1.0)]
	set avtoptilt3 [expr [lindex $av_toptilt 3]/($av_tilt_i*1.0)]
	set avtoptilt4 [expr [lindex $av_toptilt 4]/($av_tilt_i*1.0)]
	set avtoptilt5 [expr [lindex $av_toptilt 5]/($av_tilt_i*1.0)]
	set avtoptilt6 [expr [lindex $av_toptilt 6]/($av_tilt_i*1.0)]
	set avtoptilt7 [expr [lindex $av_toptilt 7]/($av_tilt_i*1.0)]
	set avtoptilt8 [expr [lindex $av_toptilt 8]/($av_tilt_i*1.0)]
	set avtoptilt9 [expr [lindex $av_toptilt 9]/($av_tilt_i*1.0)]
	set avtoptilt10 [expr [lindex $av_toptilt 10]/($av_tilt_i*1.0)] 
	set avtoptilt11 [expr [lindex $av_toptilt 11]/($av_tilt_i*1.0)]
	set avtoptilt12 [expr [lindex $av_toptilt 12]/($av_tilt_i*1.0)]
	set avtoptilt13 [expr [lindex $av_toptilt 13]/($av_tilt_i*1.0)]
	set avtoptilt14 [expr [lindex $av_toptilt 14]/($av_tilt_i*1.0)]
	set avtoptilt15 [expr [lindex $av_toptilt 15]/($av_tilt_i*1.0)]
	set avtoptilt16 [expr [lindex $av_toptilt 16]/($av_tilt_i*1.0)]
	set avtoptilt17 [expr [lindex $av_toptilt 17]/($av_tilt_i*1.0)]
	set avtoptilt18 [expr [lindex $av_toptilt 18]/($av_tilt_i*1.0)]
	set avtoptilt19 [expr [lindex $av_toptilt 19]/($av_tilt_i*1.0)]
	
	
	puts $f_tvstilt "$time $avtoptilt0 $avtoptilt1 $avtoptilt2  $avtoptilt3 $avtoptilt4 $avtoptilt5 $avtoptilt6 $avtoptilt7 $avtoptilt8 $avtoptilt9 $avtoptilt10 $avtoptilt11 $avtoptilt12  $avtoptilt13 $avtoptilt14 $avtoptilt15 $avtoptilt16 $avtoptilt17 $avtoptilt18 $avtoptilt19"
	
    } else {
	::mmsg::warn [namespace current] "can't print average tilt angle"
	flush stdout
    }
    flush $f_tvstilt
}

proc ::mbtools::analysis::setup_tiltangle { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable verbose

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_tiltangle verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_tiltangle$suffix "
    if { [file exists "$outputdir/time_vs_tiltangle$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvstilt [open "$outputdir/time_vs_tiltangle$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvstilt "\# Time angle1 angle2 angle3 angle4 angle5 angle6 angle7 angle8 angle9 angle10 angle11 angle12 angle13 angle14 angle15 angle16 angle17 angle18 angle19 angle20"
    }
}





proc ::mbtools::analysis::analyze_tiltangle {} {

    global ::mbtools::analysis::topology
   
    set molnblisttoptemp 0
    unset molnblisttoptemp


    set molnblistbottomtemp 0
    unset molnblistbottomtemp

     
#get the middle bead ids of the protein as the centre of the radial distance  

    set upperbead [::system_generation::get_upperbead]
 
    set lowerbead [::system_generation::get_lowerbead]


#turn the bead ids into molecule ids
    
    set molidupper [part $upperbead print mol]
    
    set molidlower [part $lowerbead print mol]
 
  
#get the lipid orientations
    set l_orients [lindex [analyze get_lipid_orients] 1]

    set localorients [lindex [analyze lipid_orient_order all] 1]


#a loop over the number of radial distances ($k)

###########################################
# Top layer
###########################################


    for { set k 1 } { $k < 11 } { incr k } {

	set orienttop 0

	set rcatchtop [expr (2*1.5) + $k]

	set molnblisttop [analyze nbhood planar 1 1 0 $upperbead $rcatchtop] 
	set tiltangletoptot 0

	set molnblistupper 0
	unset molnblistupper	 

	
			
#for the first layer, insert all the mol ids into the molnblistupper#
#for the second, ...layers, only the ids which are not in the previous layer are included in the molnblistupper, this is done by comparing the whole ids in the molnblisttoptemp with the new ones in the molnblisttop 

	if {$k == 1} {

	    foreach nbmol $molnblisttop {

		set thismol [part $nbmol print mol]

		if {$thismol !=  $molidupper  && [lindex $l_orients $thismol] == $orienttop} {

		lappend molnblisttoptemp $thismol

		lappend molnblistupper $thismol

	    }

	    }

	} else {
    
	    foreach nbmol $molnblisttop {

		set thismol [part $nbmol print mol]	

		if {$thismol != $molidupper && [lsearch $molnblisttoptemp $thismol] == -1 && [lindex $l_orients $thismol] == $orienttop} {
	    
		    lappend molnblistupper $thismol
		    
		    lappend molnblisttoptemp $thismol
	    
		}

	    }  
	       
	}
	

#get the distance between the protein and one of the lipid within a radial distance#
	
	set first [lindex $molnblistupper end]
	set tb [lindex [lindex $topology $first] end]
	set tbpos [bond_length $upperbead $tb]
       
	lappend disup $tbpos


#get the vectors of the molecule#

	foreach nbmol $molnblistupper {

	    set localorienttop [lindex $localorients $nbmol]

#calculate the cos theta between the bilayer vector r and the bilayer normal

	    set x [lindex $localorienttop 0]
	    set y [lindex $localorienttop 1]
	    set z [lindex $localorienttop 2]

	    set zmagnitudetop [expr { sqrt($z*$z) }]

	    set rmagnitudetop [expr { sqrt($x*$x + $y*$y + $z*$z) }]

	    set tiltangletop [expr $zmagnitudetop/$rmagnitudetop]

	    set tiltangletoptot [expr $tiltangletoptot + $tiltangletop]
    
	}

#average the angle 
    
	set tiltangletopav [expr $tiltangletoptot/[llength $molnblistupper]]

	lappend tiltangletoplist $tiltangletopav
	
    }


###############################################
#    Bottom layer
###############################################

    for { set k 1 } { $k < 11 } { incr k } {

	set orientbottom 1

	set rcatchbottom [expr (2*1.5) + $k]

	set molnblistbottom [analyze nbhood planar 1 1 0 $lowerbead $rcatchbottom] 
	
	set tiltanglebottomtot 0

	set molnblistlower 0
	unset molnblistlower


#for the first layer, insert all the mol ids into the molnblistupper#
#for the second, ...layers, only the ids which are not in the previous layer are included in the molnblistupper, this is done by comparing the whole ids in the molnblisttoptemp with the new ones in the molnblisttop 

	if {$k == 1} {

	    foreach nbmol $molnblistbottom {

		set thismol [part $nbmol print mol]

		if {$thismol !=  $molidlower  && [lindex $l_orients $thismol] == $orientbottom} {


		    lappend molnblistbottomtemp $thismol

		    lappend molnblistlower $thismol

		}
	   
	    }

	} else {
      
	    foreach nbmol $molnblistbottom {

		set thismol [part $nbmol print mol]	

		if {$thismol != $molidlower && [lsearch $molnblistbottomtemp $thismol] == -1 && [lindex $l_orients $thismol] == $orientbottom} {

		    lappend molnblistlower $thismol
		    
		    lappend molnblistbottomtemp $thismol
	    
		}

	    }
    
	}


#get the distance between the protein and a lipid within a radial distance

	set first [lindex $molnblistlower end]
	set tb [lindex [lindex $topology $first] end]
	set tbpos [bond_length $lowerbead $tb]
       
	lappend dislow $tbpos


#get the vectors of the molecule#

	foreach nbmol $molnblistlower {

	    set localorientbottom [lindex $localorients $nbmol]

#calculate the cos theta between the bilayer vector r and the bilayer normal

	    set x [lindex $localorientbottom 0]
	    set y [lindex $localorientbottom 1]
	    set z [lindex $localorientbottom 2]

	    set zmagnitudebottom [expr { sqrt($z*$z) }]

	    set rmagnitudebottom [expr { sqrt($x*$x + $y*$y + $z*$z) }]

	    set tiltanglebottom [expr $zmagnitudebottom/$rmagnitudebottom]

	    set tiltanglebottomtot [expr $tiltanglebottomtot + $tiltanglebottom]
    
	}

#average the angle 
    
	set tiltanglebottomav [expr $tiltanglebottomtot/[llength $molnblistlower]]

	lappend tiltangletoplist $tiltanglebottomav
	
    }



    return $tiltangletoplist

}

proc ::mbtools::analysis::analyze_tiltangle { printflag } {


    ::mmsg::send [namespace current] "analyzing the tilt angle of the lipids surrounding the protein"
    variable av_tilt_i
    variable av_toptilt

    set tiltangletoplist [calc_tiltangle]



    for { set k 0 } { $k < 20 } { incr k } {

	set dumy [lindex $tiltangletoplist $k]

	lappend inst_toptilt $dumy 

    }



    lset av_toptilt 0 [expr [lindex $av_toptilt 0] + [lindex $inst_toptilt 0] ]
    lset av_toptilt 1 [expr [lindex $av_toptilt 1] + [lindex $inst_toptilt 1] ]
    lset av_toptilt 2 [expr [lindex $av_toptilt 2] + [lindex $inst_toptilt 2] ]
    lset av_toptilt 3 [expr [lindex $av_toptilt 3] + [lindex $inst_toptilt 3] ]
    lset av_toptilt 4 [expr [lindex $av_toptilt 4] + [lindex $inst_toptilt 4] ]
    lset av_toptilt 5 [expr [lindex $av_toptilt 5] + [lindex $inst_toptilt 5] ]
    lset av_toptilt 6 [expr [lindex $av_toptilt 6] + [lindex $inst_toptilt 6] ]
    lset av_toptilt 7 [expr [lindex $av_toptilt 7] + [lindex $inst_toptilt 7] ]
    lset av_toptilt 8 [expr [lindex $av_toptilt 8] + [lindex $inst_toptilt 8] ]
    lset av_toptilt 9 [expr [lindex $av_toptilt 9] + [lindex $inst_toptilt 9] ]
 lset av_toptilt 10 [expr [lindex $av_toptilt 10] + [lindex $inst_toptilt 10] ]
    lset av_toptilt 11 [expr [lindex $av_toptilt 11] + [lindex $inst_toptilt 11] ]
    lset av_toptilt 12 [expr [lindex $av_toptilt 12] + [lindex $inst_toptilt 12] ]
    lset av_toptilt 13 [expr [lindex $av_toptilt 13] + [lindex $inst_toptilt 13] ]
    lset av_toptilt 14 [expr [lindex $av_toptilt 14] + [lindex $inst_toptilt 14] ]
    lset av_toptilt 15 [expr [lindex $av_toptilt 15] + [lindex $inst_toptilt 15] ]
    lset av_toptilt 16 [expr [lindex $av_toptilt 16] + [lindex $inst_toptilt 16] ]
    lset av_toptilt 17 [expr [lindex $av_toptilt 17] + [lindex $inst_toptilt 17] ]
    lset av_toptilt 18 [expr [lindex $av_toptilt 18] + [lindex $inst_toptilt 18] ]
    lset av_toptilt 19 [expr [lindex $av_toptilt 19] + [lindex $inst_toptilt 19] ]

    incr av_tilt_i

     
    if { $printflag } {

	set avtoptilt0 [expr [lindex $av_toptilt 0]/($av_tilt_i*1.0)]
	set avtoptilt1 [expr [lindex $av_toptilt 1]/($av_tilt_i*1.0)]
	set avtoptilt2 [expr [lindex $av_toptilt 2]/($av_tilt_i*1.0)]
	set avtoptilt3 [expr [lindex $av_toptilt 3]/($av_tilt_i*1.0)]
	set avtoptilt4 [expr [lindex $av_toptilt 4]/($av_tilt_i*1.0)]
	set avtoptilt5 [expr [lindex $av_toptilt 5]/($av_tilt_i*1.0)]
	set avtoptilt6 [expr [lindex $av_toptilt 6]/($av_tilt_i*1.0)]
	set avtoptilt7 [expr [lindex $av_toptilt 7]/($av_tilt_i*1.0)]
	set avtoptilt8 [expr [lindex $av_toptilt 8]/($av_tilt_i*1.0)]
	set avtoptilt9 [expr [lindex $av_toptilt 9]/($av_tilt_i*1.0)]

	set avtoptilt10 [expr [lindex $av_toptilt 10]/($av_tilt_i*1.0)]
	set avtoptilt11 [expr [lindex $av_toptilt 11]/($av_tilt_i*1.0)]
	set avtoptilt12 [expr [lindex $av_toptilt 12]/($av_tilt_i*1.0)]
	set avtoptilt13 [expr [lindex $av_toptilt 13]/($av_tilt_i*1.0)]
	set avtoptilt14 [expr [lindex $av_toptilt 14]/($av_tilt_i*1.0)]
	set avtoptilt15 [expr [lindex $av_toptilt 15]/($av_tilt_i*1.0)]
	set avtoptilt16 [expr [lindex $av_toptilt 16]/($av_tilt_i*1.0)]
	set avtoptilt17 [expr [lindex $av_toptilt 17]/($av_tilt_i*1.0)]
	set avtoptilt18 [expr [lindex $av_toptilt 18]/($av_tilt_i*1.0)]
	set avtoptilt19 [expr [lindex $av_toptilt 19]/($av_tilt_i*1.0)]


	::mmsg::send [namespace current]  "L: [lindex $inst_toptilt 0] [lindex $inst_toptilt 1] [lindex $inst_toptilt 2] [lindex $inst_toptilt 3] [lindex $inst_toptilt 4] [lindex $inst_toptilt 5] [lindex $inst_toptilt 6] [lindex $inst_toptilt 7] [lindex $inst_toptilt 8] [lindex $inst_toptilt 9] [lindex $inst_toptilt 10] [lindex $inst_toptilt 11] [lindex $inst_toptilt 12] [lindex $inst_toptilt 13] [lindex $inst_toptilt 14] [lindex $inst_toptilt 15] [lindex $inst_toptilt 16] [lindex $inst_toptilt 17] [lindex $inst_toptilt 18] [lindex $inst_toptilt 19] :: <L> $avtoptilt0 $avtoptilt1 $avtoptilt2 $avtoptilt3 $avtoptilt4 $avtoptilt5 $avtoptilt6 $avtoptilt7 $avtoptilt8  $avtoptilt9 $avtoptilt10 $avtoptilt11 $avtoptilt12 $avtoptilt13 $avtoptilt14 $avtoptilt15 $avtoptilt16 $avtoptilt17 $avtoptilt18 $avtoptilt19 "

	flush stdout
    }	

  
    ::mmsg::debug [namespace current] "done"

}