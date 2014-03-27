# Copyright (C) 2012,2013 The ESPResSo project
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

proc oif_init {} {	
    # define oif variables
        
    global oif_n_objects   
    set oif_n_objects 0
    # variable that denotes the number of objects in the simulation

    global oif_n_templates
    set oif_n_templates 0
    # variable that denotes the number of templates of objects in the simulation

    global list oif_nparticles
    # list of the numbers of nodes for each oif object

    global list oif_ntriangles
    # list of the numbers of triangles for each oif object

    global list oif_nedges
    # list with the numbers of edges for each oif object

    global oif_first_bond_id 
    set oif_first_bond_id 0
    # variable that denotes the ID of bond that can be created next

    global oif_first_part_id
    set oif_first_part_id 0
    # variable that denotes the ID of particle that can be created next

    global oif_first_inter_id 
    set oif_first_inter_id 0
    # denotes the ID of interaction, that can be used for the next oif_template
	
    global list oif_templates
    # 2D list of existing templates of oif objects containing rows with data. One row describes parameters of one object template. One row consists of X elements. The order of the elements is crucial. Structure of one row: num_of_particles, num_of_edges, num_of_triangles, ks, start_id_of_ks_interaction, kb, start_id_of_kb_interaction, kal, start_id_of_kal_interaction, kag, start_id_of_kag_interaction, kv, start_id_of_kv_interaction

    global list oif_objects
    # 2D list of existing objects containing rows with data. One row describes parameters of one object. One row consists of X elements. The order of the elements is crucial. Structure of one row: template_id, part_type, part_mass

    global list oif_template_particles
    # list containing triplets of coordinates of all particles of all templates that will be used when creating objects. First num_of_particles of them belong to the first template. Next num_of_particles belong to the second template, etc. Note: these are not actual espresso particles, but rather node coordinates that will serve as template for the actual objects. 

    global list oif_template_edges
    # list containing pairs of particle IDs. The ordering works the same way as for particles. First num_of_edges edges belong to the first template, next num_of_edges belong to the second template, etc.

    global list oif_template_triangles
    # list containing triplets of particle IDs. Again, in this list, triangles of all templates are arranged consequently and oif_template element num_of_triangles is used to determine which triangles belong to which template.
	
    global list oif_template_bending_incidences
    # a list of bending incidences (4 particle IDs in each row; number of rows equals number of rows in oif_template_edges). These are used so that during the creation of bending bonds one does not have to recalculate the bending incidences repeatedly. Number of bending bonds of a given template is equal to number of edges of this template, so oif_template element num_of_edges is used to determine which incidences belong to which object.

    global list oif_object_starting_particles
    # a list of particle indices denoting where each object's particles start in the whole list of particles

    global list oif_template_starting_triangles
    # a list of triangle indices denoting where each object's triangles start in the whole list of triangles

}

proc oif_info { } {
	global oif_n_objects
        global oif_n_templates
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_first_bond_id
	global oif_first_part_id
	global oif_first_inter_id
	global oif_templates
	global oif_template_particles
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_object_starting_particles
	global oif_template_starting_triangles

	puts " "
	puts "*************************************"
	puts "*                                   *"
	puts "*       Info about oif objects      *"
	puts "*                                   *"
	puts "*************************************"
	
	puts "oif_n_objects: $oif_n_objects"
        puts "oif_n_templates: $oif_n_templates"

	puts "oif_nparticles"
	foreach item $oif_nparticles {
		puts "$item"	
	}
	puts "oif_ntriangles"
	foreach item $oif_ntriangles {
		puts "$item"	
	}
	puts "oif_nedges"
	foreach item $oif_nedges {
		puts "$item"	
	}
	
	puts "oif_first_bond_id: $oif_first_bond_id"
	puts "oif_first_part_id: $oif_first_part_id"
        puts "oif_first_inter_id: $oif_first_inter_id"

	puts "oif_templates"
	foreach item $oif_templates {
		puts "$item"	
	}

        puts "oif_template_particles"
	foreach item $oif_template_particles {
		puts "$item"	
	}

        puts "oif_template_edges"
	foreach item $oif_template_edges {
		puts "$item"	
	}

        puts "oif_template_triangles"
	foreach item $oif_template_triangles {
		puts "$item"	
	}

        puts "oif_template_bending_incidences"
	foreach item $oif_template_bending_incidences {
		puts "$item"	
	}
	
        puts "oif_object_starting_particles"
	foreach item $oif_object_starting_particles {
		puts "$item"	
	}

        puts " "
	puts " "
	puts " "

}

proc get_n_triangle {a b c n} {
	upvar $a ga
	upvar $b gb 
	upvar $c gc 
	upvar $n gn 
	# now ga refers to the "global" variable with name $a
	set P1x [lindex $ga 0]
	set P1y [lindex $ga 1]
	set P1z [lindex $ga 2]
	set P2x [lindex $gb 0]
	set P2y [lindex $gb 1]
	set P2z [lindex $gb 2]
	set P3x [lindex $gc 0]
	set P3y [lindex $gc 1]
	set P3z [lindex $gc 2]

	set nx [expr ($P1y - $P2y)*($P3z - $P2z) - ($P1z - $P2z)*($P3y - $P2y)]
	set ny [expr ($P1z - $P2z)*($P3x - $P2x) - ($P1x - $P2x)*($P3z - $P2z)]
	set nz [expr ($P1x - $P2x)*($P3y - $P2y) - ($P1y - $P2y)*($P3x - $P2x)]

	set gn [list $nx $ny $nz]
}

proc area_triangle {a b c} {
	# Computes the area of triangle P1,P2,P3 by computing the crossproduct P1P2 x P1P3 and taking the half of its norm
	upvar $a ga
	upvar $b gb 
	upvar $c gc 

	set area 0
	set n {0 0 0} 
	
	get_n_triangle ga gb gc n
	set nx [lindex $n 0]
	set ny [lindex $n 1]
	set nz [lindex $n 2]

	set area [expr 0.5*sqrt($nx*$nx + $ny*$ny + $nz*$nz)]
	return $area
}

proc angle_btw_triangles {P1 P2 P3 P4 phi} {
	upvar $P1 gP1
	upvar $P2 gP2
	upvar $P3 gP3
	upvar $P4 gP4
	upvar $phi gphi
	set n1 {-1 -1 -1}
	set n2 {-1 -1 -1}
	set pi 3.1415926535897931
	
	get_n_triangle gP1 gP2 gP3 n1
	get_n_triangle gP3 gP2 gP4 n2
	
	set n1x [lindex $n1 0]
	set n1y [lindex $n1 1]
	set n1z [lindex $n1 2]
	set n2x [lindex $n2 0]
	set n2y [lindex $n2 1]
	set n2z [lindex $n2 2]
	
	set tmp11 [expr $n1x*$n2x + $n1y*$n2y + $n1z*$n2z]
	set tmp11 [expr $tmp11*abs($tmp11)]
	set tmp22 [expr $n1x*$n1x + $n1y*$n1y + $n1z*$n1z]
	set tmp33 [expr $n2x*$n2x + $n2y*$n2y + $n2z*$n2z]
	set tmp11 [expr $tmp11/($tmp22*$tmp33)] 

	if { $tmp11 > 0 } {
		set tmp11 [expr sqrt($tmp11)]
		} else {
		set tmp11 [expr - sqrt(- $tmp11)]
	}		

	if {$tmp11 >= 1.} { set tmp11 0.0 } elseif { $tmp11 <= -1.} { set tmp11 $pi }
	
	set gphi [expr $pi - acos($tmp11)]

	set P1x [lindex $gP1 0]
	set P1y [lindex $gP1 1]
	set P1z [lindex $gP1 2]
	set P4x [lindex $gP4 0]
	set P4y [lindex $gP4 1]
	set P4z [lindex $gP4 2]
	
	set tmp11 [expr -($n1x*$P1x + $n1y*$P1y + $n1z*$P1z)]
	if { [expr $n1x*$P4x + $n1y*$P4y + $n1z*$P4z + $tmp11] < 0 } { set gphi [expr 2*$pi - $gphi] }
}

proc discard_epsilon {x} {
	set res 0.0
	if { $x > -0.0000000001 &&  $x < 0.0000000001 } { 
		set res 0.0 
	} else { 
		set res $x 
	}
	return $res
}

proc oif_create_template { args } {

# access global variables defined in oif_init
	global oif_first_inter_id
	global oif_n_templates	
	global oif_templates
	global oif_template_particles
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_template_starting_triangles

#--------------------------------------------------------------------------------------------
# count the number of arguments and initialize
        set n_args 0
	foreach arg $args {
		incr n_args
	}
	if { $n_args == 0 } {
		puts "Mandatory arguments are nodes file, triangles file"
		return 0
	}

	set stretch_X 1.0
	set stretch_Y 1.0
	set stretch_Z 1.0
	set filenamenodes ""
	set filenametriangles ""
	set ks 0.0
	set kb 0.0
	set kal 0.0
	set kag 0.0
	set kv 0.0
	set check_output 0
	set template_id -1

#--------------------------------------------------------------------------------------------
# reading the arguments. some of them are mandatory. we check for the mandatory arguments at the end of this section
	set pos 0
	while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"nodes-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set filenamenodes [lindex $args $pos]
				incr pos
			}
			"ks" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set ks [lindex $args $pos]
				incr pos
			}
			"kb" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set kb [lindex $args $pos]
				incr pos
			}
			"kal" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set kal [lindex $args $pos]
				incr pos
			}
			"kag" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set kag [lindex $args $pos]
				incr pos
			}
			"kv" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set kv [lindex $args $pos]
				incr pos
			}
			"triangles-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set filenametriangles [lindex $args $pos]
				incr pos
			}
			"stretch" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set stretch_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set stretch_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set stretch_Z [lindex $args $pos]
				incr pos
			}
		        "check" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set check_output [lindex $args $pos]
				incr pos
			}
		        "template-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set template_id [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error" 
				set pos $n_args
			}
		}  
	}

#--------------------------------------------------------------------------------------------
# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $filenamenodes == "" } { set mandatory 0 }
	if { $filenametriangles == "" } { set mandatory 0 }
	if { $template_id == -1 } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for template creator"  
		return
	}
	
#--------------------------------------------------------------------------------------------	
# checking whether correct template-id was given
	if { $template_id != $oif_n_templates } {
		puts "error: Cannot create a template with template_id $template_id,"
		puts "	because either one such template already exists or the template_id does not follow consecutively."
		puts "	The next available template_id is $oif_n_templates."
		exit	
	}

#--------------------------------------------------------------------------------------------
# set files for output check
	if {$check_output == 1} {

		set bondS "TMP/noGbondsStretching$template_id"
		set bondB "TMP/noGbondsBending$template_id"
		set bondAlocal "TMP/noGbondsAreaLocal$template_id"
		set bondAglobal "TMP/noGbondsAreaGlobal$template_id"
		set bondV "TMP/noGbondsVolume$template_id" 
	}

#--------------------------------------------------------------------------------------------
# check: output all parameters:
	puts "The following oif_template will be created:"
	puts "	nodes-file: 	$filenamenodes"
	puts "	triangles-file: $filenametriangles"
	puts "	ks:		$ks"
	puts "	kb:		$kb"
	puts "	kal: 		$kal"
	puts "	kag: 		$kag"
	puts "	kv:		$kv"
	puts "	stretch: 	$stretch_X $stretch_Y $stretch_Z"
#--------------------------------------------------------------------------------------------
#reading nodes
	# read the number of lines in nodes files
	set fp [open $filenamenodes r]
	set file_data [read $fp]
	close $fp
	set data [split $file_data "\n"]

	# template must be stretched first
	set mesh_nnodes 0

	foreach line $data {
		if { [llength $line] == 3 } {
			set mesh_nodes($mesh_nnodes,0) [expr $stretch_X*[lindex $line 0]]
			set mesh_nodes($mesh_nnodes,1) [expr $stretch_Y*[lindex $line 1]]
			set mesh_nodes($mesh_nnodes,2) [expr $stretch_Z*[lindex $line 2]]
			
			# mesh_nodes is a 2D-array with three coordinates for each node (each node is one line) 
		        # node $X coordinate y is accessed by $mesh_nodes($X,1)
			incr mesh_nnodes
		        
	        list temp_particle
	        lappend temp_particle [expr $stretch_X*[lindex $line 0]]
			lappend temp_particle [expr $stretch_Y*[lindex $line 1]]
			lappend temp_particle [expr $stretch_Z*[lindex $line 2]]
	        lappend oif_template_particles $temp_particle
	        # update global variable
	        set temp_particle [lreplace $temp_particle 0 2]
	        # clear list
		}	    
	}

#--------------------------------------------------------------------------------------------
#reading triangles
	# read the number of lines in triangle file
	set fp [open $filenametriangles r]
	set file_data [read $fp]
	close $fp
	set data [split $file_data "\n"]

	set mesh_ntriangles 0
	# setting the first triangle id for the currently created template
	if {$template_id > 0 } { 
		lappend oif_template_starting_triangles [llength $oif_template_triangles] 
	} else { lappend oif_template_starting_triangles 0 }

	foreach line $data {
		if { [llength $line] == 3 } {
			set mesh_triangles($mesh_ntriangles,0) [lindex $line 0]
			set mesh_triangles($mesh_ntriangles,1) [lindex $line 1]
			set mesh_triangles($mesh_ntriangles,2) [lindex $line 2]
			incr mesh_ntriangles
		        
	        list temp_triangle
	        lappend temp_triangle [lindex $line 0]
			lappend temp_triangle [lindex $line 1]
			lappend temp_triangle [lindex $line 2]
	        lappend oif_template_triangles $temp_triangle
	        # update global variable
	        set temp_triangle [lreplace $temp_triangle 0 2]
	        # clear list
		}
	}
	
#--------------------------------------------------------------------------------------------
# check data extracted from input files:
	puts "Data extracted from the input files:"
	puts "	nnodes: 	$mesh_nnodes"
	puts "	ntriangles: 	$mesh_ntriangles"

#--------------------------------------------------------------------------------------------
# creating the list of edges	

	set mesh_nedges 0
	
	for {set i 0} {$i < $mesh_ntriangles} {incr i} {
	    # take a triangle and copy the nodes of the triangle to pa,pb,pc (point A, point B, point C)
	    set pa $mesh_triangles($i,0)
	    set pb $mesh_triangles($i,1)
	    set pc $mesh_triangles($i,2)
	    set is 0	

	    for {set j 0} {$j < $mesh_nedges} {incr j} {
		# Check if the edge AB or BA is in the current list of edges
		if {$mesh_edges($j,0) == $pa && $mesh_edges($j,1) == $pb} { set is 1}
		if {$mesh_edges($j,1) == $pa && $mesh_edges($j,0) == $pb} { set is 1}   
		}

	    if { $is == 0} {  
		# If AB nor BA is in the list then add the edge AB to the list
		set mesh_edges($mesh_nedges,0) $pa
		set mesh_edges($mesh_nedges,1) $pb
		incr mesh_nedges
		}

	    set is 0

	    for {set j 0} {$j < $mesh_nedges} {incr j} {
		# Check if the edge BC or CB is in the current list of edges
		if {$mesh_edges($j,0) == $pb && $mesh_edges($j,1) == $pc} { set is 1}
		if {$mesh_edges($j,1) == $pb && $mesh_edges($j,0) == $pc} { set is 1}  
		}

	    if {$is == 0} {  
		# If BC nor CB is in the list then add the edge BC to the list
		set mesh_edges($mesh_nedges,0) $pb
		set mesh_edges($mesh_nedges,1) $pc
		incr mesh_nedges
		}

	    set is 0

	    for {set j 0} {$j < $mesh_nedges} {incr j} {
		# Check if the edge AC or CA is in the current list of edges
		if {$mesh_edges($j,0) == $pa && $mesh_edges($j,1) == $pc} { set is 1}
		if {$mesh_edges($j,1) == $pa && $mesh_edges($j,0) == $pc} { set is 1}
		}

	    if {$is == 0} {     
		# If AC nor CA is in the list then add the edge AC to the list
		set mesh_edges($mesh_nedges,0) $pa
		set mesh_edges($mesh_nedges,1) $pc
		incr mesh_nedges
		}
	}

	for {set i 0} {$i < $mesh_nedges} {incr i} {
	    list temp_edge
		        lappend temp_edge $mesh_edges($i,0)
			lappend temp_edge $mesh_edges($i,1)
		        lappend oif_template_edges $temp_edge
		        # update global variable
		        set temp_edge [lreplace $temp_edge 0 2]
		        # clear list
	}

#----------------------------------------------------------------------------------------
# prepare template data for global variable oif_template 
	list template
	lappend template $mesh_nnodes
	lappend template $mesh_nedges
	lappend template $mesh_ntriangles
#-----------------------------------------------------------------------------------------
# generating stretching force interactions 
	lappend template $ks
	if { $ks == 0.0} {
	    set start_id_of_ks_interaction -1 
	} else {
	    set start_id_of_ks_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_ks_interaction

	if { $ks != 0.0} {

	    if {$check_output == 1} { 
	       	set fpart [open $bondS "w"]
	    }

	    puts "generating stretching force interactions"
	    # Stretching is coupled to the edges   
	    set n_StrInter $mesh_nedges
	    set oif_first_inter_id [expr $oif_first_inter_id+$n_StrInter]

	    set dist 0
	    for {set i 0} {$i < $n_StrInter} {incr i} {
		for {set k 0} {$k < 3} {incr k} {		
		    set p1($k) [expr $mesh_nodes($mesh_edges($i,0),$k)]
		    set p2($k) [expr $mesh_nodes($mesh_edges($i,1),$k)]
		}
		# We need to compute the distance btw the vertices
		set dist [expr sqrt(($p1(0)-$p2(0))*($p1(0)-$p2(0)) + ($p1(1)-$p2(1))*($p1(1)-$p2(1)) + ($p1(2)-$p2(2))*($p1(2)-$p2(2)))] 
		inter [expr $start_id_of_ks_interaction + $i] stretching_force [format %e $dist] [format %e $ks]
		if {$check_output == 1} { 
		    puts $fpart "inter [expr $start_id_of_ks_interaction + $i] stretching_force [format %e $dist] [format %e $ks]"
		}
	    }
	    if {$check_output == 1} { 
		close $fpart	   
	    }
	}	

#-----------------------------------------------------------------------------------------
# generating bending force interactions 
	lappend template $kb
	if { $kb == 0.0} {
	    set start_id_of_kb_interaction -1 
	} else {
	    set start_id_of_kb_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_kb_interaction

	if { $kb != 0.0} {
	
	    if {$check_output == 1} { 
	       	set fpart [open $bondB "w"]
	    }
	    puts "generating bending force interactions"
	    set n_BenInter $mesh_nedges
	    # Bending is coupled to the angles between triangles sharing the same edge
	    set oif_first_inter_id [expr $oif_first_inter_id + $n_BenInter]
	
	    set phi 0.0
	    for { set i 0} { $i < $n_BenInter} {incr i} { 
			#Run over all edges
			set p2id $mesh_edges($i,0)
			#Put IDs of points to p2id,p3id
			set p3id $mesh_edges($i,1) 
			for { set k 0} {$k < 3} {incr k} {	
				#Put coordinates of the edges' points
				set p2($k) $mesh_nodes($p2id,$k)
				set p3($k) $mesh_nodes($p3id,$k)
			}
			
			set detected 0
				#Number of detected triangles with current edge common
				# Algorithm is as follows: we run over all triangles and check whether two vertices are those from current edge. If we find such triangle, we put the ID of the third vertex to p1id and moreover we check if the orientation p1id, p2id p3id is the same as was in the triangle list (meaning, that we found one of the following three triples in the triangle list: p1id, p2id, p3id or p2id, p3id, p1id or p3id, p1id, p2id). If we have the same orientation, we set orient = 1, otherwise orient = -1.
				# Then we go further looking for the second triangle. The second triangle should have the opposite orientation.
				# The normal of the first triangle will be P1P2 x P1P3, of the second triangle will be P2P4 x P2P3
			
			set orient 0
				# We have two different orientations. I set orientation to 1 if the first triangle is oriented p1id
	
			for { set k 0} { $k < $mesh_ntriangles} {incr k} { 
				# Run over all triangles and determine the two triangles with the common current edge
				if {$k < $mesh_ntriangles && $mesh_triangles($k,0) == $p2id && $mesh_triangles($k,1) == $p3id} {  
					if {$detected == 0} { 
							#if no triangle was detected
						set p1id $mesh_triangles($k,2) 
						set detected 1 
						set orient 1
					} else { 
							# if already one triangle was detected - then also quit the k-loop
						set p4id $mesh_triangles($k,2) 
						set k $mesh_ntriangles
					} 
				}
				if {$k < $mesh_ntriangles && $mesh_triangles($k,1) == $p2id && $mesh_triangles($k,2) == $p3id} {
					 if {$detected == 0} {
						 set p1id $mesh_triangles($k,0)
						 set detected 1
						 set orient 1
					} else {
						set p4id $mesh_triangles($k,0)
						set k $mesh_ntriangles
					}
				}
				if {$k < $mesh_ntriangles && $mesh_triangles($k,2) == $p2id && $mesh_triangles($k,0) == $p3id} {
					 if {$detected == 0} {
						 set p1id $mesh_triangles($k,1)
						 set detected 1 
						 set orient 1
					} else {
						set p4id $mesh_triangles($k,1)
						set k $mesh_ntriangles
					} 
				}
				if {$k < $mesh_ntriangles && $mesh_triangles($k,1) == $p2id && $mesh_triangles($k,0) == $p3id} {
					 if {$detected == 0} {
						 set p1id $mesh_triangles($k,2)
						 set detected 1 
						 set orient -1
					 } else {
						set p4id $mesh_triangles($k,2)
						set k $mesh_ntriangles
					}
				}
				if {$k < $mesh_ntriangles && $mesh_triangles($k,2) == $p2id && $mesh_triangles($k,1) == $p3id} {
					 if {$detected == 0} {
						 set p1id $mesh_triangles($k,0)
						 set detected 1 
						 set orient -1
					} else {
						set p4id $mesh_triangles($k,0)
						set k $mesh_ntriangles
					}
				}
				if {$k < $mesh_ntriangles && $mesh_triangles($k,0) == $p2id && $mesh_triangles($k,2) == $p3id} {
					 if {$detected == 0} {
						 set p1id $mesh_triangles($k,1)
						 set detected 1 
						 set orient -1
					} else {
						set p4id $mesh_triangles($k,1)
						set k $mesh_ntriangles
					}
				}
			}
			if {$orient == 1} {
				set tmp22 $p1id
				set p1id $p4id
				set p4id $tmp22
			} 
				#This is to have the correct orientation
			list P1
			lappend P1 $mesh_nodes($p1id,0)
			lappend P1 $mesh_nodes($p1id,1)
			lappend P1 $mesh_nodes($p1id,2)
			list P2
			lappend P2 $mesh_nodes($p2id,0)
			lappend P2 $mesh_nodes($p2id,1)
			lappend P2 $mesh_nodes($p2id,2)
			list P3
			lappend P3 $mesh_nodes($p3id,0)
			lappend P3 $mesh_nodes($p3id,1)
			lappend P3 $mesh_nodes($p3id,2)
			list P4
			lappend P4 $mesh_nodes($p4id,0)
			lappend P4 $mesh_nodes($p4id,1)
			lappend P4 $mesh_nodes($p4id,2)
	
			list temp_incidence
		        lappend temp_incidence $p1id
			lappend temp_incidence $p2id
			lappend temp_incidence $p3id
			lappend temp_incidence $p4id
		        lappend oif_template_bending_incidences $temp_incidence
		        # update global variable
		        set temp_incidence [lreplace $temp_incidence 0 3] 
		        # clear list

			angle_btw_triangles P1 P2 P3 P4 phi

			# to replace lists with empty lists so they do not grow
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
			set P4 [lreplace $P4 0 2]

			inter [expr $start_id_of_kb_interaction + $i] bending_force [format %e $phi] [format %e $kb]
			if {$check_output == 1} { 
			    puts $fpart "inter [expr $start_id_of_kb_interaction + $i] bending_force [format %e $phi] [format %e $kb]"
			}
		}
	
	    if {$check_output == 1} { 
		close $fpart	   
	    }
	}
#-----------------------------------------------------------------------------------------------
# generation of local area force interactions
	lappend template $kal
	if { $kal == 0.0} {
	    set start_id_of_kal_interaction -1 
	} else {
	    set start_id_of_kal_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_kal_interaction

	if {$kal != 0.0} {

	    if {$check_output == 1} { 
	       	set fpart [open $bondAlocal "w"]
	    }

	    puts "generating local area force interactions"
	    set n_localAreaInter $mesh_ntriangles
	    # Area is coupled to the triangles
	    set oif_first_inter_id [expr $oif_first_inter_id + $n_localAreaInter]
	
	    set area 0.0
	    for {set i 0} {$i < $n_localAreaInter} {incr i} {
		for {set k 0} {$k < 3} {incr k} {
		    list P1
		    lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
		    list P2
		    lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
		    list P3
		    lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
		}
		set area [area_triangle P1 P2 P3]
		# to replace lists with empty lists so they do not grow
		set P1 [lreplace $P1 0 2]
		set P2 [lreplace $P2 0 2]
		set P3 [lreplace $P3 0 2]
			
		inter [expr $start_id_of_kal_interaction + $i] area_force_local $area $kal
		if {$check_output == 1} { 
		    puts $fpart "inter [expr $start_id_of_kal_interaction + $i] area_force_local $area $kal"
		}
	    }
	    if {$check_output == 1} { 
		close $fpart	   
	    }
	}
#--------------------------------------------------------------------------------------------
# generation of global area force interactions
	lappend template $kag
	if { $kag == 0.0} {
	    set start_id_of_kag_interaction -1 
	} else {
	    set start_id_of_kag_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_kag_interaction

	if {$kag != 0.0} {
	    
	    if {$check_output == 1} { 
	       	set fpart [open $bondAglobal "w"]
	    }
	    
	    puts "generating global area force interactions"
	    set n_globalAreaInter 1
	    set oif_first_inter_id [expr $oif_first_inter_id + $n_globalAreaInter]
	
	    set area 0.0
	    set gl_area 0.0

	    for {set i 0} {$i < $mesh_ntriangles} {incr i} {
		for {set k 0} {$k < 3} {incr k} {		
		    list P1
		    lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
		    list P2
		    lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
		    list P3
		    lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
		}
		set area [area_triangle P1 P2 P3]
		# to replace lists with empty lists so they do not grow
		set P1 [lreplace $P1 0 2]
		set P2 [lreplace $P2 0 2]
		set P3 [lreplace $P3 0 2]
		set gl_area [expr $gl_area + $area]
	    }
	    inter $start_id_of_kag_interaction area_force_global $gl_area $kag
	    if {$check_output == 1} { 
		puts $fpart "inter $start_id_of_kag_interaction area_force_global $gl_area $kag"
		close $fpart
	    }
	}
#--------------------------------------------------------------------------------------------	
# generation of volume force interactions
	lappend template $kv
	if { $kv == 0.0} {
	    set start_id_of_kv_interaction -1 
	} else {
	    set start_id_of_kv_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_kv_interaction

	if {$kv != 0.0} {

	    if {$check_output == 1} { 
	       	set fpart [open $bondV "w"]
	    }
	    
	    puts "generating volume force interactions"
	    set n_VolumeInter 1
	    set oif_first_inter_id [expr $oif_first_inter_id + $n_VolumeInter]
	
	    set area 0.0
	    set volume 0.0 
	    set hz 0.0
	    list norm
	    set dn 0.0
	    set drmax 0.0
	    set drtemp 0.0
		
	    for {set i 0} {$i < $mesh_ntriangles} {incr i} {
		for {set k 0} {$k < 3} {incr k} {		
		    list P1
		    lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
		    list P2
		    lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
		    list P3
		    lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
		}
		set P10 [lindex $P1 0]
		set P11 [lindex $P1 1]
		set P12 [lindex $P1 2]
		set P20 [lindex $P2 0]
		set P21 [lindex $P2 1]
		set P22 [lindex $P2 2]
		set P30 [lindex $P3 0]
		set P31 [lindex $P3 1]
		set P32 [lindex $P3 2]
		set drtemp [expr sqrt(($P10 - $P20)*($P10 - $P20) + ($P11 - $P21)*($P11 - $P21) + ($P12 - $P22)*($P12 - $P22))]
		# distance P1P2
		if {$drmax < $drtemp} { set drmax $drtemp }
		set drtemp [expr sqrt(($P10 - $P30)*($P10 - $P30) + ($P11 - $P31)*($P11 - $P31) + ($P12 - $P32)*($P12 - $P32))]
		# distance P1P3
		if {$drmax < $drtemp} { set drmax $drtemp }
			
		set area [area_triangle P1 P2 P3]
		get_n_triangle P1 P2 P3 norm
		set norm0 [lindex $norm 0]
		set norm1 [lindex $norm 1]
		set norm2 [lindex $norm 2]
		set dn [expr sqrt($norm0*$norm0 + $norm1*$norm1 + $norm2*$norm2)]
		set hz [expr 1.0/3.0*($P12 + $P22 + $P32)]
		set volume [expr $volume + $area*$norm2/$dn*$hz]

		# to replace lists with empty lists so they do not grow
		set P1 [lreplace $P1 0 2]
	       	set P2 [lreplace $P2 0 2]
		set P3 [lreplace $P3 0 2]
	    }
	    inter $start_id_of_kv_interaction volume_force $volume $kv
	    if {$check_output == 1} { 
		puts $fpart "inter $start_id_of_kv_interaction volume_force $volume $kv"
		close $fpart
	    }
	}

#--------------------------------------------------------------------------------------------
# update global oif-variables        
	lappend oif_templates $template
	incr oif_n_templates
} 

#--------------------------------------------------------------------------------------------
proc oif_add_object { args } {
	
# access global variables defined in init_objects_in_fluid
	global oif_n_objects
	global oif_n_templates
	global oif_first_bond_id
	global oif_first_part_id
	global oif_first_inter_id
	global oif_templates
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_template_particles
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_object_starting_particles
	global oif_objects

#--------------------------------------------------------------------------------------------
# counts the number of arguments and initialize
	set n_args 0
	foreach arg $args {
		incr n_args
	}
	if { $n_args == 0 } {
		puts "Mandatory arguments are object template, origin, particle type, particle mol"
		return 0
	}

	set rotate_X 0
	set rotate_Y 0
	set rotate_Z 0
	set origin_X 0
	set origin_Y 0
	set origin_Z 0
	set part_type -1
	set part_mol -1
	set part_mass 1.0
	set template_id -1
	set object_id -1
	set check_output 0

#--------------------------------------------------------------------------------------------
# reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
	set pos 0
	while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"mass" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set part_mass [lindex $args $pos]
				incr pos
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set object_id [lindex $args $pos]
				set part_mol $object_id
				incr pos
			}
			"part-type" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set part_type [lindex $args $pos]
				incr pos
			}
			"template-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set template_id [lindex $args $pos]
				incr pos
			}
			"rotate" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error rot 1"
					break
				}
				set rotate_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error rot 2"
					break
				}
				set rotate_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error rot 3"
					break
				}
				set rotate_Z [lindex $args $pos]
				incr pos
			}
			"origin" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error origin 1"
					break
				}
				set origin_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error  origin 1"
					break
				}
				set origin_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error  origin 1"
					break
				}
				set origin_Z [lindex $args $pos]
				incr pos
			}
			"check" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error check"
					break
				}
				set check_output [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error default" 
				set pos $n_args
			}
		}  
	}


#--------------------------------------------------------------------------------------------
# checking wheter all mandatory arguments have been given
	set mandatory 1
	if { $origin_X == 0 &&  $origin_Y == 0 &&  $origin_Z == 0 } { set mandatory 0 }
	if { $part_type == "-1" } { set mandatory 0 }
	if { $template_id == "-1" } { set mandatory 0 }
	if { $object_id == "-1" } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for creating object" 
		return
	}
	
#--------------------------------------------------------------------------------------------	
# checking whether correct object-id has been given
	if { $object_id != $oif_n_objects } {
		puts "error: Cannot create an object with object_id $object_id,"
		puts "	because either one such object already exists or the object_id does not follow consecutively."
		puts "	The next available object_id is $oif_n_objects."
		exit	
	}
	
#--------------------------------------------------------------------------------------------
# set files for output check
	if {$check_output == 1} {
		set createPart "TMP/noGcreatePart$object_id"
		set partS "TMP/noGpartStretching$object_id"
		set partB "TMP/noGpartBending$object_id"
		set partAlocal "TMP/noGpartAreaLocal$object_id"
		set partAglobal "TMP/noGpartAreaGlobal$object_id"
		set partV "TMP/noGpartVolume$object_id"
	}

#--------------------------------------------------------------------------------------------
# Check: output all parameters:
	puts "The following oif_object will be created:"
	puts "	template_id:	$template_id"
	puts "	mass:	 	$part_mass"
	puts "	particle_type: 	$part_type"
	puts "	mol: 		$part_mol"
	puts "	rotate: 	$rotate_X $rotate_Y $rotate_Z"
	puts "	origin: 	$origin_X $origin_Y $origin_Z"

#----------------------------------------------------------------------------------------
# prepare object data for global variable oif_objects
	list object_data
	lappend object_data $template_id
	lappend object_data $part_type
	lappend object_data $part_mass
#--------------------------------------------------------------------------------------------
# find and read the template for the given object

 	# there are 13 elements in one row of oif_templates:
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	# get the data form oif_templates list
	set template [lindex $oif_templates $template_id]
	
	set mesh_nnodes [lindex $template 0]
	set mesh_nedges [lindex $template 1]
	set mesh_ntriangles [lindex $template 2]	

	set start_id_of_ks_interactions [lindex $template 4]
	set start_id_of_kb_interactions [lindex $template 6]
	set start_id_of_kal_interactions [lindex $template 8]
	set start_id_of_kag_interactions [lindex $template 10]
	set start_id_of_kv_interactions [lindex $template 12]

	set start_id_of_particles 0
	set start_id_of_edges 0
	set start_id_of_triangles 0
	for {set i 0} {$i < $template_id} {incr i} {
	    set start_id_of_particles [expr $start_id_of_particles + [lindex [lindex $oif_templates $i] 0]]
	    set start_id_of_edges [expr $start_id_of_edges + [lindex [lindex $oif_templates $i] 1]]
	    set start_id_of_triangles [expr $start_id_of_triangles + [lindex [lindex $oif_templates $i] 2]]
	}
	set start_id_of_bending_incidences $start_id_of_edges
	 
 	# recover particles of the given template  
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
	    set node_triplet [lindex $oif_template_particles [expr $start_id_of_particles+$i]]
	    set mesh_nodes($i,0) [lindex $node_triplet 0] 
	    set mesh_nodes($i,1) [lindex $node_triplet 1]
	    set mesh_nodes($i,2) [lindex $node_triplet 2]	    
	}

	# recover triangles of the given template

	for {set i 0} {$i < $mesh_ntriangles} {incr i} {
	    set triangle_triplet [lindex $oif_template_triangles [expr $start_id_of_triangles+$i]]
	    set mesh_triangles($i,0) [lindex $triangle_triplet 0]
	    set mesh_triangles($i,1) [lindex $triangle_triplet 1]
	    set mesh_triangles($i,2) [lindex $triangle_triplet 2]
	}

	# recover edges of the given template
	for {set i 0} {$i < $mesh_nedges} {incr i} {
	    set edge_pair [lindex $oif_template_edges [expr $start_id_of_edges+$i]]
	    set mesh_edges($i,0) [lindex $edge_pair 0]
	    set mesh_edges($i,1) [lindex $edge_pair 1]
	}

	# recover bending incidences of the given template	
	# template data:
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set kb_from_template [lindex $template 5]
	if { $kb_from_template != 0.0 } {
		for {set i 0} {$i < $mesh_nedges} {incr i} {
		    set bending_quartet [lindex $oif_template_bending_incidences [expr $start_id_of_bending_incidences+$i]]
		    set bending_incidences($i,0) [lindex $bending_quartet 0]
		    set bending_incidences($i,1) [lindex $bending_quartet 1]
		    set bending_incidences($i,2) [lindex $bending_quartet 2]
		    set bending_incidences($i,3) [lindex $bending_quartet 3]
		}
	}
#--------------------------------------------------------------------------------------------
# some variables for rotation
	set ca [expr cos($rotate_X)];
	set sa [expr sin($rotate_X)];
	set cb [expr cos($rotate_Y)];
	set sb [expr sin($rotate_Y)];
	set cc [expr cos($rotate_Z)];
	set sc [expr sin($rotate_Z)];
	set rotation(0,0) [expr $cb*$cc]
	set rotation(0,1) [expr $sa * $sb * $cc - $ca * $sc]
	set rotation(0,2) [expr $sc * $sa + $cc * $sb * $ca]
	set rotation(1,0) [expr $cb * $sc]
	set rotation(1,1) [expr $ca * $cc - $sa * $sb * $sc ]
	set rotation(1,2) [expr $sc * $sb * $ca - $cc * $sa]
	set rotation(2,0) [expr -$sb]
	set rotation(2,1) [expr $cb * $sa ]
	set rotation(2,2) [expr $ca * $cb  ]
	
# rotation of nodes:
 	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set xx [discard_epsilon [expr $rotation(0,0)*$mesh_nodes($i,0) + $rotation(0,1)*$mesh_nodes($i,1) + $rotation(0,2)*$mesh_nodes($i,2)]]
		set yy [discard_epsilon [expr $rotation(1,0)*$mesh_nodes($i,0) + $rotation(1,1)*$mesh_nodes($i,1) + $rotation(1,2)*$mesh_nodes($i,2)]]
		set zz [discard_epsilon [expr $rotation(2,0)*$mesh_nodes($i,0) + $rotation(2,1)*$mesh_nodes($i,1) + $rotation(2,2)*$mesh_nodes($i,2)]]

		set mesh_nodes($i,0) $xx;
		set mesh_nodes($i,1) $yy;
		set mesh_nodes($i,2) $zz;
	}
#--------------------------------------------------------------------------------------------
# shift to origin:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set mesh_nodes($i,0) [expr $mesh_nodes($i,0) + $origin_X]
		set mesh_nodes($i,1) [expr $mesh_nodes($i,1) + $origin_Y]
		set mesh_nodes($i,2) [expr $mesh_nodes($i,2) + $origin_Z]
	}

#--------------------------------------------------------------------------------------------  	
#
#
#         ESPRESSO object creation
#	
# 	
#--------------------------------------------------------------------------------------------
# generating particles:
	puts "generating particles"
	if {$check_output == 1} { set f [open $createPart "w"] }
	set i $oif_first_part_id
	
	# remember where this object's particles start:
	lappend oif_object_starting_particles $oif_first_part_id

	for {set i $oif_first_part_id} {$i < [expr $mesh_nnodes + $oif_first_part_id]} {incr i} {
		part $i pos [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],0)]] [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],1)]] [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],2)]] type $part_type mol $part_mol mass $part_mass
		if {$check_output == 1} { puts $f [format "part $i pos %e %e %e type $part_type mol $part_mol mass $part_mass" [expr $mesh_nodes([expr $i - $oif_first_part_id],0)] [expr $mesh_nodes([expr $i - $oif_first_part_id],1)] [expr $mesh_nodes([expr $i - $oif_first_part_id],2)]] }
		}  
	set oif_first_part_id [expr $oif_first_part_id + $mesh_nnodes]
	if {$check_output == 1} { close $f }

#--------------------------------------------------------------------------------------------
# generation of stretching force bonds:
	# template data
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set ks_from_template [lindex $template 3]
	if { $ks_from_template != 0.0 } {
		if {$check_output == 1} { 
		       	set fpart [open $partS "w"]
		}
	   	puts "generating stretching force bonds"
		set firstID_StrInter $start_id_of_ks_interactions 
		# Stretching is coupled to the edges   
		set n_StrBond $mesh_nedges 
		set oif_first_bond_id [expr $oif_first_bond_id+$n_StrBond]
		  
	   	for {set i 0} {$i < $n_StrBond} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
		    part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrInter + $i] [expr $mesh_edges($i,1) + $firstPartId]
		    if {$check_output == 1} { 
			puts $fpart "part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrInter + $i] [expr $mesh_edges($i,1) + $firstPartId]"
		    }
	   	}
		if {$check_output == 1} { 
		   close $fpart	   
	    }
	}

#--------------------------------------------------------------------------------------------		
# generation of bending force bonds:
	# template data
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set kb_from_template [lindex $template 5]
	if { $kb_from_template != 0.0 } {
		if {$check_output == 1} { 
		    set fpart [open $partB "w"]
		}
			
		puts "generating bending force bonds"
		set firstID_BenInter $start_id_of_kb_interactions
		set n_BenBond $mesh_nedges
		# Bending is coupled to the angles between triangles sharing the same edge
		set oif_first_bond_id [expr $oif_first_bond_id + $n_BenBond]
		
		for {set i 0} {$i < $n_BenBond} {incr i} {
	
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
	
		    part [expr $bending_incidences($i,1) + $firstPartId] bond [expr $firstID_BenInter + $i] [expr $bending_incidences($i,0) + $firstPartId] [expr $bending_incidences($i,2) + $firstPartId] [expr $bending_incidences($i,3) + $firstPartId]
	
		    if {$check_output == 1} { 
			puts $fpart "part [expr $bending_incidences($i,1) + $firstPartId] bond [expr $firstID_BenInter + $i] [expr $bending_incidences($i,0) + $firstPartId] [expr $bending_incidences($i,2) + $firstPartId] [expr $bending_incidences($i,3) + $firstPartId]"
		    }
		}
		if {$check_output == 1} { 
		    close $fpart	
		}
	}
#--------------------------------------------------------------------------------------------
# generation of local area force bonds
	# template data
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set kal_from_template [lindex $template 7]
	if { $kal_from_template != 0.0 } {
	
		if {$check_output == 1} { 
			set fpart [open $partAlocal "w"]
		}
		puts "generating local area force bonds"
		set firstID_localAreaInter $start_id_of_kal_interactions
		set n_localAreaBond $mesh_ntriangles
		# Area is coupled to the triangles
		set oif_first_bond_id [expr $oif_first_bond_id + $n_localAreaBond]
		
	       	for {set i 0} {$i < $n_localAreaBond} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
		    part [expr $mesh_triangles($i,0) + $firstPartId] bond [expr $firstID_localAreaInter + $i] [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]
		    if {$check_output == 1} { 
			puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond [expr $firstID_localAreaInter + $i] [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]"
		    }
	       	}
		if {$check_output == 1} {
				close $fpart
		}
	}
#--------------------------------------------------------------------------------------------
# generation of global area force bonds
	# template data
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set kag_from_template [lindex $template 9]
	if { $kag_from_template != 0.0 } {
	
		if {$check_output == 1} { 
		    set fpart [open $partAglobal "w"]
	       	}
		puts "generating global area force bonds"
		set firstID_globalAreaInter $start_id_of_kag_interactions
		set n_globalAreaBond 1
		set oif_first_bond_id [expr $oif_first_bond_id + $n_globalAreaBond]
	
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
		    if {$check_output == 1} { 
			puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalAreaInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]" 
		    }
		    part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalAreaInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]
		}
		if {$check_output == 1} { 
		    close $fpart
	       	}
	}
#--------------------------------------------------------------------------------------------	
# generation of volume force bonds
	# template data
	# nnodes, nparticles, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter
	set kv_from_template [lindex $template 11]
	if { $kv_from_template != 0.0 } {
		if {$check_output == 1} { 
		    set fpart [open $partV "w"]
	       	}
		puts "generating volume force bonds"
		set firstID_VolumeInter $start_id_of_kv_interactions
		set n_VolumeBond 1
		set oif_first_bond_id [expr $oif_first_bond_id + $n_VolumeBond]
		
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
		    if {$check_output == 1} { 
			puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_VolumeInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) +$firstPartId]" 
		    }
		    part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_VolumeInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) +$firstPartId]
		}
		if {$check_output == 1} { 
		    close $fpart		
		}
	}
#--------------------------------------------------------------------------------------------	
# update global oif-variables
	incr oif_n_objects
	lappend oif_nparticles $mesh_nnodes
	lappend oif_ntriangles $mesh_ntriangles
	lappend oif_nedges $mesh_nedges
	lappend oif_objects $object_data
	
}

proc oif_object_set { args } {
	# acces global variables defined in init_objects_in_fluid.tcl:
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_triangles
	global oif_firstBondId
	global oif_firstPartId
	global oif_objects
	global oif_object_starting_particles
	global oif_template_starting_triangles

	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are: object's ID (id)"
		return 0
	}
	
	set objectID -1
	set force 0
	set mesh_nodes_file ""
	set origin 0
	set kill_motion 0
	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"force" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set forceX [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set forceY [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set forceZ [lindex $args $pos]
				incr pos
				set force 1
			}
			"origin" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set originX [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set originY [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set originZ [lindex $args $pos]
				incr pos
				set origin 1
			}
			"mesh-nodes" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"kill-motion" {  
				incr pos
				set kill_motion 1
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error s" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for print_oif_object" 
		return
	}

# setting the ext_force using part command
	if { $force == 1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			part $iii ext_force $forceX $forceY $forceZ
		}
	}

	if { $kill_motion == 1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			part $iii fix
		}
	}

	
	if { $mesh_nodes_file != ""} {
		# first find the current center of the object
		set centerX 0
		set centerY 0
		set centerZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set centerX [expr $centerX + [lindex $coords 0]]
			set centerY [expr $centerY + [lindex $coords 1]]
			set centerZ [expr $centerZ + [lindex $coords 2]]
		}
		set centerX [expr $centerX/$nnode]
		set centerY [expr $centerY/$nnode]
		set centerZ [expr $centerZ/$nnode]

		#then set new coordinates
		set fp [open $mesh_nodes_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		
		set n_lines 0
		foreach line $data {
			if { [llength $line] == 3 } {
				incr n_lines
			}
		}
		
		if { $n_lines == $nnode } {
			set iii $firstPartId
			foreach line $data {
				if { [llength $line] == 3 } {
					set posX [expr [lindex $line 0] + $centerX]
					set posY [expr [lindex $line 1] + $centerY]
					set posZ [expr [lindex $line 2] + $centerZ]
					part $iii pos $posX $posY $posZ
					incr iii
				}
			}
		} else {
			puts "Number of lines in $mesh_nodes_file is not the same as claimed number of mesh nodes for the current object."
		}
	}

	if { $origin == 1} {
		# first find the current center of the object
		set centerX 0
		set centerY 0
		set centerZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set centerX [expr $centerX + [lindex $coords 0]]
			set centerY [expr $centerY + [lindex $coords 1]]
			set centerZ [expr $centerZ + [lindex $coords 2]]
		}
		set centerX [expr $centerX/$nnode]
		set centerY [expr $centerY/$nnode]
		set centerZ [expr $centerZ/$nnode]

		#then set new coordinates
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set newposX [expr $originX - $centerX + [lindex $coords 0]]
			set newposY [expr $originY - $centerY + [lindex $coords 1]]
			set newposZ [expr $originZ - $centerZ + [lindex $coords 2]]
			part $iii pos $newposX $newposY $newposZ
		}
	}


	
}


proc oif_object_output { args } {
	# acces global variables defined in init_objects_in_fluid.tcl:
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_triangles
	global oif_object_starting_particles
	global oif_firstBondId
	global oif_firstPartId
	global oif_firstTriangleId
	global oif_objects
	global oif_template_starting_triangles
	global oif_template_triangles
	
	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are: object's ID (id)"
		return 0
	}
	
	set objectID -1
	set vtk_pos_file ""
	set vtk_aff_file ""
	set mesh_nodes_file ""

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"vtk-pos" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set vtk_pos_file [lindex $args $pos]
				incr pos
			}
			"mesh-nodes" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error rets" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for oif_object_print" 
		return
	}

	if { $vtk_pos_file != ""} {
		set part [open $vtk_pos_file "w"]
		puts $part "# vtk DataFile Version 3.0"
		puts $part "Data"
		puts $part "ASCII"
		puts $part "DATASET POLYDATA"
		puts $part "POINTS [lindex $oif_nparticles $objectID] float"
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			puts $part "[part $iii print pos]"
		}
		puts $part "TRIANGLE_STRIPS [lindex $oif_ntriangles $objectID] [expr 4*[lindex $oif_ntriangles $objectID]]"
		# gets info about the current object
		set object_data [lindex $oif_objects $objectID] 
		# extracts the id of the object's template
		set object_template_id [lindex $object_data 0]
		# extracts the starting position of the triangles for the template_id
		set firstTriangleId [lindex $oif_template_starting_triangles $object_template_id]
		for { set iii $firstTriangleId } { $iii < [expr $firstTriangleId + [lindex $oif_ntriangles $objectID]] } { incr iii } {
			puts $part "3 [lindex $oif_template_triangles $iii]"
		}
		close $part
	}

	if { $mesh_nodes_file != ""} {
		set centerX 0
		set centerY 0
		set centerZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set centerX [expr $centerX + [lindex $coords 0]]
			set centerY [expr $centerY + [lindex $coords 1]]
			set centerZ [expr $centerZ + [lindex $coords 2]]
		}
		set centerX [expr $centerX/$nnode]
		set centerY [expr $centerY/$nnode]
		set centerZ [expr $centerZ/$nnode]

		set part [open $mesh_nodes_file "w"]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set posX [expr - $centerX + [lindex $coords 0]]
			set posY [expr - $centerY + [lindex $coords 1]]
			set posZ [expr - $centerZ + [lindex $coords 2]]
			puts $part "$posX $posY $posZ"
		}
		close $part
		set coords [list $centerX $centerY $centerZ]
		return $coords
	}

}


proc oif_object_analyze { args } {
	# acces global variables defined in init_objects_in_fluid.tcl:
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_triangles
	global oif_object_starting_particles
	global oif_firstBondId
	global oif_firstPartId
	global oif_firstTriangleId
	global oif_objects
	global oif_template_starting_triangles
	global oif_template_triangles

	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are: object's ID (id)"
		return 0
	}
	
	set objectID -1
	set center_location 0
	set pos_bounds ""
	set volume -1
	set surface -1
	set energy ""
	set velocity 0
	set affinity ""

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"origin" {  
				incr pos
				set center_location 1
			}
			"velocity" {  
				incr pos
				set velocity 1
			}
			"volume" {  
				incr pos
				set volume 0.0
			}
			"surface-area" {  
				incr pos
				set surface 0.0
			}
			"pos-bounds" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set pos_bounds [lindex $args $pos]
				incr pos
			}
			"energy" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set energy [lindex $args $pos]
				incr pos
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in reading the arguments of oif_object_analyze" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for oif_object_print" 
		return
	}

	if { $center_location == 1 } {
		set centerX 0
		set centerY 0
		set centerZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set centerX [expr $centerX + [lindex $coords 0]]
			set centerY [expr $centerY + [lindex $coords 1]]
			set centerZ [expr $centerZ + [lindex $coords 2]]
		}
		set centerX [expr $centerX/$nnode]
		set centerY [expr $centerY/$nnode]
		set centerZ [expr $centerZ/$nnode]
		set coords [list $centerX $centerY $centerZ]
		return $coords
	}

	if { $velocity == 1 } {
		set velX 0
		set velY 0
		set velZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set vel [part $iii print v]
			set velX [expr $velX + [lindex $vel 0]]
			set velY [expr $velY + [lindex $vel 1]]
			set velZ [expr $velZ + [lindex $vel 2]]
		}
		set velX [expr $velX/$nnode]
		set velY [expr $velY/$nnode]
		set velZ [expr $velZ/$nnode]
		set vel [list $velX $velY $velZ]
		return $vel
	}



	if {$volume != -1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		set ntriangles [lindex $oif_ntriangles $objectID]
		set counter 0
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set mesh_nodes($counter,0) [lindex $coords 0]
			set mesh_nodes($counter,1) [lindex $coords 1]
			set mesh_nodes($counter,2) [lindex $coords 2]
			incr counter
		}
		set counter 0
		# gets info about the current object
		set object_data [lindex $oif_objects $objectID] 
		# extracts the id of the object's template
		set object_template_id [lindex $object_data 0]
		# extracts the starting position of the triangles for the template_id
		set firstTriangleId [lindex $oif_template_starting_triangles $object_template_id]
		for { set iii $firstTriangleId } { $iii < [expr $firstTriangleId + $ntriangles] } { incr iii } {
			set incidence_line [lindex $oif_template_triangles $iii]
			set mesh_triangles($counter,0) [lindex $incidence_line 0]
			set mesh_triangles($counter,1) [lindex $incidence_line 1]
			set mesh_triangles($counter,2) [lindex $incidence_line 2]
			incr counter
		}
		set area 0.0
		set volume 0.0 
		set hz 0.0
		list norm
		set dn 0.0
		set drmax 0.0
		set drtemp 0.0
		
		for {set i 0} {$i < $ntriangles} {incr i} {
			for {set k 0} {$k < 3} {incr k} {		
				list P1
				lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
				list P2
				lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
				list P3
				lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
			}
			set P10 [lindex $P1 0]
			set P11 [lindex $P1 1]
			set P12 [lindex $P1 2]
			set P20 [lindex $P2 0]
			set P21 [lindex $P2 1]
			set P22 [lindex $P2 2]
			set P30 [lindex $P3 0]
			set P31 [lindex $P3 1]
			set P32 [lindex $P3 2]
			set drtemp [expr sqrt(($P10 - $P20)*($P10 - $P20) + ($P11 - $P21)*($P11 - $P21) + ($P12 - $P22)*($P12 - $P22))]
				# distance P1P2
			if {$drmax < $drtemp} { set drmax $drtemp }
			set drtemp [expr sqrt(($P10 - $P30)*($P10 - $P30) + ($P11 - $P31)*($P11 - $P31) + ($P12 - $P32)*($P12 - $P32))]
				# distance P1P3
			if {$drmax < $drtemp} { set drmax $drtemp }
			
			set area [area_triangle P1 P2 P3]
			get_n_triangle P1 P2 P3 norm
			set norm0 [lindex $norm 0]
			set norm1 [lindex $norm 1]
			set norm2 [lindex $norm 2]
			set dn [expr sqrt($norm0*$norm0 + $norm1*$norm1 + $norm2*$norm2)]
			set hz [expr 1.0/3.0*($P12 + $P22 + $P32)]
			set volume [expr $volume + $area*$norm2/$dn*$hz]
			# to replace lists with empty lists so they do not grow
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
		}
		return $volume
	}
	
	if {$surface != -1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		set ntriangles [lindex $oif_ntriangles $objectID]
		set counter 0
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set mesh_nodes($counter,0) [lindex $coords 0]
			set mesh_nodes($counter,1) [lindex $coords 1]
			set mesh_nodes($counter,2) [lindex $coords 2]
			incr counter
		}
		set counter 0
		# gets info about the current object
		set object_data [lindex $oif_objects $objectID] 
		# extracts the id of the object's template
		set object_template_id [lindex $object_data 0]
		# extracts the starting position of the triangles for the template_id
		set firstTriangleId [lindex $oif_template_starting_triangles $object_template_id]
		for { set iii $firstTriangleId } { $iii < [expr $firstTriangleId + $ntriangles] } { incr iii } {
			set incidence_line [lindex $oif_template_triangles $iii]
			set mesh_triangles($counter,0) [lindex $incidence_line 0]
			set mesh_triangles($counter,1) [lindex $incidence_line 1]
			set mesh_triangles($counter,2) [lindex $incidence_line 2]
			incr counter
		}
		set area 0.0
		for {set i 0} {$i < $ntriangles} {incr i} {
			for {set k 0} {$k < 3} {incr k} {		
				list P1
				lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
				list P2
				lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
				list P3
				lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
			}
			set curr_area [area_triangle P1 P2 P3]
			set area [expr $area + $curr_area]
			# to replace lists with empty lists so they do not grow
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
		}
		return $area
	}
	
	if { $pos_bounds != ""} {
		set Zmin 10000000
		set Zmax -10000000
		set Xmin 10000000
		set Xmax -10000000
		set Ymin 10000000
		set Ymax -10000000

		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set coords [part $iii print pos]
			set Xcoord [lindex $coords 0]
			set Ycoord [lindex $coords 1]
			set Zcoord [lindex $coords 2]
			if { $Xcoord < $Xmin } { set Xmin $Xcoord }
			if { $Xcoord > $Xmax } { set Xmax $Xcoord }
			if { $Ycoord < $Ymin } { set Ymin $Ycoord }
			if { $Ycoord > $Ymax } { set Ymax $Ycoord }
			if { $Zcoord < $Zmin } { set Zmin $Zcoord }
			if { $Zcoord > $Zmax } { set Zmax $Zcoord }
		}
		if { $pos_bounds == "all" } { 
			set result [list $Xmin $Xmax $Ymin $Ymax $Zmin $Zmax]
			return $result
		}
		if { $pos_bounds == "z-min" } { set result $Zmin }
		if { $pos_bounds == "z-max" } { set result $Zmax }
		if { $pos_bounds == "x-min" } { set result $Xmin }
		if { $pos_bounds == "x-max" } { set result $Xmax }
		if { $pos_bounds == "y-min" } { set result $Ymin }
		if { $pos_bounds == "y-max" } { set result $Ymax }

		return $result
	}

	if { $energy != ""} {
		# TO BE IMPLEMENTED ......
		set result -1
		if { $energy == "ks" } {
			set result 0
		}
		if { $energy == "kb" } {
			set result 0
		}
		if { $energy == "kal" } {
			set result 0
		}
		if { $energy == "kag" } {
			set result 0
		}
		if { $energy == "kv" } {
			set result 0
		}
		if { $energy == "aff" } {
			set result 0
		}
		if { $result == -1 } {
			puts "Argument of energy computation must be one of these: ks, kb, kal, kag, kv, aff"
		}
		return $result
	}

   
}


proc oif_mesh_analyze { args } {
	# acces global variables defined in init_objects_in_fluid.tcl:
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_triangles
	global oif_object_starting_particles
	global oif_firstBondId
	global oif_firstPartId

	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are: nodes-file, triangles-file"
		return 0
	}
	
	set mesh_nodes_file ""
	set mesh_triangles_file ""
	set orientation 0
	set corr_file ""
	set method -1
	set flip -1
	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"nodes-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"triangles-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error 4"
					break
				}
				set mesh_triangles_file [lindex $args $pos]
				incr pos
			}
			"orientation" {
				incr pos
				set orientation 1
			}
			"repair" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set corr_file [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set method [lindex $args $pos]
				incr pos
			}
			"flip" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set corr_file [lindex $args $pos]
				set flip 1
				incr pos
			}
			default { 
				puts "error s" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $mesh_nodes_file == "" } { set mandatory 0 }
	if { $mesh_triangles_file == "" } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with mandatory arguments for oif_check_mesh" 
		return
	}

# checking the orientation 
	if { $orientation == 1} {
		# Method 1: First method runs through all triangles and checks whether point (0,0,0) is on the same side of all planes given by the triangles. This check is not applicable for non-convex objects.
		set fp [open $mesh_nodes_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		set mesh_nnodes 0
		foreach line $data {
			if { [llength $line] == 3 } {
				set mesh_nodes($mesh_nnodes,0) [lindex $line 0]
				set mesh_nodes($mesh_nnodes,1) [lindex $line 1]
				set mesh_nodes($mesh_nnodes,2) [lindex $line 2]
					# mesh_nodes is an 2D-array with  three coordinates for each node. (each node is one line) you can access node $X coordinate y  by $mesh_nodes($X,1)
				incr mesh_nnodes
			}
		}
		
		set fp [open $mesh_triangles_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		set mesh_ntriangles 0
		foreach line $data {
			if { [llength $line] == 3 } {
				set mesh_triangles($mesh_ntriangles,0) [lindex $line 0]
				set mesh_triangles($mesh_ntriangles,1) [lindex $line 1]
				set mesh_triangles($mesh_ntriangles,2) [lindex $line 2]
				lappend oif_triangles $line
				incr mesh_ntriangles
			}
		}
		
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			for {set k 0} {$k < 3} {incr k} {		
				list P1
				lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
				list P2
				lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
				list P3
				lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
			}
			set P10 [lindex $P1 0]
			set P11 [lindex $P1 1]
			set P12 [lindex $P1 2]
			#set P20 [lindex $P2 0]
			#set P21 [lindex $P2 1]
			#set P22 [lindex $P2 2]
			#set P30 [lindex $P3 0]
			#set P31 [lindex $P3 1]
			#set P32 [lindex $P3 2]
			#set drtemp [expr sqrt(($P10 - $P20)*($P10 - $P20) + ($P11 - $P21)*($P11 - $P21) + ($P12 - $P22)*($P12 - $P22))]
				## distance P1P2
			#if {$drmax < $drtemp} { set drmax $drtemp }
			#set drtemp [expr sqrt(($P10 - $P30)*($P10 - $P30) + ($P11 - $P31)*($P11 - $P31) + ($P12 - $P32)*($P12 - $P32))]
				## distance P1P3
			#if {$drmax < $drtemp} { set drmax $drtemp }
			
			#set area [area_triangle P1 P2 P3]
	
			get_n_triangle P1 P2 P3 norm
			set norm0 [lindex $norm 0]
			set norm1 [lindex $norm 1]
			set norm2 [lindex $norm 2]
			
			# first compute "d" from the normal equation of the triangle plane
			set tmp_d [expr -($P10*$norm0 + $P11*$norm1 + $P12*$norm2)]
			# when coordinates of the origin are placed into the normal equation of the triangle plane, you should always get a positive number.
			#set tmp_res [expr $origin_X*$norm0 + $origin_Y*$norm1 + $origin_Z*$norm2 + $tmp_d]
#			puts "tmp_d: $tmp_d"
			if { $tmp_d >= 0 } {
				puts "something is wrong with the orientation of the triangles!!! at line $i"
			}
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]

		}
		
		# Method 2: Second check controls whether all couples of triangles with the same edge have the same orientation. Not implementedn yet

	}
	
	
	if { $method == 1} {
		# Method 1: First method runs through all triangles and checks whether point (0,0,0) is on the same side of all planes given by the triangles. This check is not applicable for non-convex objects.
		set fp [open $mesh_nodes_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		set mesh_nnodes 0
		foreach line $data {
			if { [llength $line] == 3 } {
				set mesh_nodes($mesh_nnodes,0) [lindex $line 0]
				set mesh_nodes($mesh_nnodes,1) [lindex $line 1]
				set mesh_nodes($mesh_nnodes,2) [lindex $line 2]
					# mesh_nodes is an 2D-array with  three coordinates for each node. (each node is one line) you can access node $X coordinate y  by $mesh_nodes($X,1)
				incr mesh_nnodes
			}
		}
		
		set fp [open $mesh_triangles_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		set mesh_ntriangles 0
		foreach line $data {
			if { [llength $line] == 3 } {
				set mesh_triangles($mesh_ntriangles,0) [lindex $line 0]
				set mesh_triangles($mesh_ntriangles,1) [lindex $line 1]
				set mesh_triangles($mesh_ntriangles,2) [lindex $line 2]
				lappend oif_triangles $line
				incr mesh_ntriangles
			}
		}
		
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			for {set k 0} {$k < 3} {incr k} {		
				list P1
				lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
				list P2
				lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
				list P3
				lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
			}
			set P10 [lindex $P1 0]
			set P11 [lindex $P1 1]
			set P12 [lindex $P1 2]
			get_n_triangle P1 P2 P3 norm
			set norm0 [lindex $norm 0]
			set norm1 [lindex $norm 1]
			set norm2 [lindex $norm 2]
			
			# first compute "d" from the normal equation of the triangle plane
			set tmp_d [expr -($P10*$norm0 + $P11*$norm1 + $P12*$norm2)]
			# when coordinates of the origin are placed into the normal equation of the triangle plane, you should always get a positive number.
			#set tmp_res [expr $origin_X*$norm0 + $origin_Y*$norm1 + $origin_Z*$norm2 + $tmp_d]
#			puts "tmp_d: $tmp_d"
			if { $tmp_d >= 0 } {
				set tmp_id $mesh_triangles($i,1) 
				set mesh_triangles($i,1) $mesh_triangles($i,2) 
				set mesh_triangles($i,2) $tmp_id
			}
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
		}

		set fp [open $corr_file w]
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			puts $fp "$mesh_triangles($i,0) $mesh_triangles($i,1) $mesh_triangles($i,2)"
		}		
		close $fp
	}


	if { $method == 2} {
		# Method 2: Second check controls whether all couples of triangles with the same edge have the same orientation. Not implementedn yet
		# TO BE IMPLEMENTED
	}

	if { $flip == 1} {
		set fcorr [open $corr_file w]
		set fp [open $mesh_triangles_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		foreach line $data {
			if { [llength $line] == 3 } {
				set v0 [lindex $line 0]
				set v1 [lindex $line 2]
				set v2 [lindex $line 1]
				puts $fcorr "$v0 $v1 $v2"
			}
		}
		close $fcorr
	}

	
}


