# Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
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

source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "CONSTRAINTS"
require_feature "PARTIAL_PERIODIC"
require_feature "ADRESS" off

if { [catch {

    puts "----------------------------------------"
    puts "- Testcase maggs_varying_epsilon.tcl running on [format %02d [setmd n_nodes]] nodes: -"
    puts "----------------------------------------"

	set accepted_error 1e-2

    # Tuning parameters
    set time_step   0.01
    set skin        0.3


    # Interaction parameters
    #############################################################

    # Coulomb
    set bjerrum       5.0
    set f_mass        0.01
    set mesh          24

    set tcl_precision 7

    #############################################################
    #  Setup System                                             #
    #############################################################

    # setup new configuration
    set box_l 10.0
    setmd box_l $box_l $box_l $box_l

	set posx [expr $box_l / 2.0 ]
	set posy1 [expr $box_l / 4.0 ]
	set posy2 [expr $box_l * 3.0/4.0 ]
	set posz1 [expr $box_l / 3.0 ]
	set posz2 [expr $box_l * 2.0 / 3.0 ]

	part 0 pos $posx $posy1 $posz1 v 0.0 0.0 0.0 f 0.0 0.0 0.0 q -1.0
	part 1 pos $posx $posy2 $posz2 v 0.0 0.0 0.0 f 0.0 0.0 0.0 q 1.0

#    setmd max_num_cells $max_cells3d

    setmd time_step $time_step
    setmd skin      $skin      
    thermostat off
    integrate 0

    # memd requires domain decompostion with no verlet lists
    cellsystem domain_decomposition -no_verlet_list
    inter coulomb $bjerrum memd $f_mass $mesh

	# set permittivity to 80 at      1/4 * box_l  <  z  <  3/4 * box_l
	for {set x 0} {$x < $mesh} {incr x} {
		for {set y 0} {$y < $mesh} {incr y} {
			for {set z [expr int(floor($mesh/4.0))] } {$z < [expr int(floor(3.0*$mesh/4.0))] } {incr z} {
				inter coulomb $bjerrum memd localeps node $x $y $z dir X eps 80.0
				inter coulomb $bjerrum memd localeps node $x $y $z dir Y eps 80.0
				inter coulomb $bjerrum memd localeps node $x $y $z dir Z eps 80.0
			}
		}
	}


    #############################################################
    #      Integration                                          #
    #############################################################


	integrate 50
	set forcez [ expr [ lindex [ part 0 print f ] 2 ] / 5.413 ]

	# compare test results with ELC result

    # unset MEMD
    inter coulomb 0.0
    # set cell size
    set box_l_z 5.0
#    setmd box_l $box_l $box_l $box_l_z
    # set periodicity for MMM2D
#    setmd periodic 1 1 0
    # MMM2D requires layered cell system
#	cellsystem layered 3 -no_verlet_list

	# reset particle coordinates, speeds and forces
	part 0 pos $posx $posy1 $posz1 v 0.0 0.0 0.0 f 0.0 0.0 0.0 q -1.0
	part 1 pos $posx $posy2 $posz2 v 0.0 0.0 0.0 f 0.0 0.0 0.0 q 1.0

	# set up ELC-MMM2D
#	inter coulomb $bjerrum mmm2d 1e-5 dielectric 80.0 1.0 80.0

	integrate 50
	set forcez_elc [ lindex [ part 0 print f ] 2 ]


    puts "\n\nTestcase maggs_varying_epsilon.tcl finished.\n"
    if { [ expr abs($forcez - $forcez_elc) ] > $accepted_error } {
    	puts "MEMD result: $forcez"
    	puts "ELC result: $forcez_elc"
		error "MEMD: deviation is too large.\n"
    } else {
		puts "MEMD with varying epsilon seems to have no errors.\n"
    } 

} res ] } {
    error_exit $res
}

exit 0