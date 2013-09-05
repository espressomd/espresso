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

proc init_objects_in_fluid {} {
	
	# define oif-variables
	global oif_n_objects
	set oif_n_objects 0
		# denotes the number of oif objects
	global list oif_nnode
		# list with the numbers of nodes for each oif object
	global list oif_ntriangles
		# list with the numbers of triangles for each oif object
	global list oif_nedges
		# list with the numbers of edges for each oif object
	global oif_firstBondId 
	set oif_firstBondId 0
	    # denotes the ID of bonds, that can be used for the next oif 
	global oif_firstPartId
	set oif_firstPartId 0
	    # denotes the ID of particle, that can be used for the next oif 
}

proc oif_info { } {
	global oif_n_objects
	global oif_nnode
	global oif_ntriangles
	global oif_nedges
	global oif_firstBondId
	global oif_firstPartId

	puts " "
	puts "*************************************"
	puts "*                                   *"
	puts "*       Info about oif objects      *"
	puts "*                                   *"
	puts "*************************************"
	
	puts "oif_n_objects: $oif_n_objects"
	puts "oif_nnode"
	foreach item $oif_nnode {
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
	
	puts "oif_firstBondId: $oif_firstBondId"
	puts "oif_firstPartId: $oif_firstPartId"

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

	set area 0;
	set n {0 0 0} 
	
	get_n_triangle ga gb gc n;
	set nx [lindex $n 0]
	set ny [lindex $n 1]
	set nz [lindex $n 2]

	set area [expr 0.5*sqrt($nx*$nx + $ny*$ny + $nz*$nz)]
	return $area;
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
	
	set tmp11 [expr $n1x*$n2x + $n1y*$n2y + $n1z*$n2z];
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

proc add_oif_object { args } {
	# acces global variables defined in init_objects_in_fluid.tcl:
	global oif_n_objects
	global oif_nnode
	global oif_ntriangles
	global oif_nedges
	global oif_firstBondId
	global oif_firstPartId

	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are origin, nodes file, triangles file, particle type, particle mol"
		return 0
	}
	set rotate_X 0
	set rotate_Y 0
	set rotate_Z 0
	set origin_X 0
	set origin_Y 0
	set origin_Z 0
	set stretch_X 1.0
	set stretch_Y 1.0
	set stretch_Z 1.0
	set filenamenodes ""
	set filenametriangles ""
	set ks 0.0
	set kslin 0.0
	set kb 0.0
	set kal 0.0
	set kag 0.0
	set kv 0.0
	set part_type -1
	set part_mol -1
	set part_mass 1
	set check_output 0

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"nodesfile" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set filenamenodes [lindex $args $pos]
				incr pos
			}
			"mass" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set part_mass [lindex $args $pos]
				incr pos
			}
			"type" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set part_type [lindex $args $pos]
				incr pos
			}
			"mol" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set part_mol [lindex $args $pos]
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
			"kslin" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set kslin [lindex $args $pos]
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
			"trianglesfile" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set filenametriangles [lindex $args $pos]
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
			"rotate" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set rotate_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set rotate_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set rotate_Z [lindex $args $pos]
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
			"origin" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set origin_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					breaken given
				}
				set origin_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error"
					break
				}
				set origin_Z [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error" 
				set pos $n_args
			}
		}  
	}


# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $origin_X == 0 &&  $origin_Y == 0 &&  $origin_Z == 0 } { set mandatory 0 }
	if { $filenamenodes == "" } { set mandatory 0 }
	if { $filenametriangles == "" } { set mandatory 0 }
	if { $part_type == "-1" } { set mandatory 0 }
	if { $part_mol == "-1" } { set mandatory 0 }
	
	if { $mandatory == 0 } { 
		puts "Something went wrong with some mandatory arguments for bond_generator" 
		return
	}

	if {$check_output == 1} {
		set createPart TMP/noGcreatePart;
		set bondS TMP/noGbondsStretching
		set bondSLIN TMP/noGbondsStretchlin
		set bondB TMP/noGbondsBending
		set bondAlocal TMP/noGbondsAreaLocal
		set bondAglobal TMP/noGbondsAreaGlobal
		set bondV TMP/noGbondsVolume
		set bondVA TMP/noGbondsVolumeAreaGlobal
		set partS TMP/noGpartStretching
		set partSLIN TMP/noGpartStretchlin
		set partB TMP/noGpartBending
		set partAlocal TMP/noGpartAreaLocal
		set partAglobal TMP/noGpartAreaGlobal
		set partV TMP/noGpartVolume
		set partVA TMP/noGpartVolumeAreaGlobal
	}


# Check: output all parameters:
	puts "The following oif-object has been created:"
	puts "	nodesfile: 		$filenamenodes"
	puts "	trianglesfile: 	$filenametriangles"
	puts "	mass:	 	$part_mass"
	puts "	type 		$part_type"
	puts "	mol: 		$part_mol"
	puts "	ks:		$ks"
	puts "	kb:		$kb"
	puts "	kal: 		$kal"
	puts "	kag: 		$kag"
	puts "	kv:		$kv"
	puts "	rotate: 	$rotate_X $rotate_Y $rotate_Z"
	puts "	stretch: 	$stretch_X $stretch_Y $stretch_Z"
	puts "	origin: 	$origin_X $origin_Y $origin_Z"


# read the number of lines in nodes files
	set fp [open $filenamenodes r]
	set file_data [read $fp]
	close $fp
	set data [split $file_data "\n"]
	set mesh_nnodes 0;
	foreach line $data {
		if { [llength $line] == 3 } {
			set mesh_nodes($mesh_nnodes,0) [lindex $line 0]
			set mesh_nodes($mesh_nnodes,1) [lindex $line 1]
			set mesh_nodes($mesh_nnodes,2) [lindex $line 2]
				# mesh_nodes is an 2D-array with  three coordinates for each node. (each node is one line) you can access node $X coordinate y  by $mesh_nodes($X,1)
			incr mesh_nnodes
		}
	}

# read the number of lines in triangle file
	set fp [open $filenametriangles r]
	set file_data [read $fp]
	close $fp

	set data [split $file_data "\n"]
	set mesh_ntriangles 0;
	foreach line $data {
		if { [llength $line] == 3 } {
			set mesh_triangles($mesh_ntriangles,0) [lindex $line 0]
			set mesh_triangles($mesh_ntriangles,1) [lindex $line 1]
			set mesh_triangles($mesh_ntriangles,2) [lindex $line 2]
			incr mesh_ntriangles
		}
	}

# Check for data extracted from input files:
	puts "Data extracted from the input files:"
	puts "	nnodes: 	$mesh_nnodes"
	puts "	ntriangles: 	$mesh_ntriangles"

# basic checks for the mesh:
	# all triangles should have correct orientation 
	# TO BE IMPLEMENTED

	# stretching of the object-in-fluid:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set mesh_nodes($i,0) [expr $mesh_nodes($i,0)*$stretch_X]
		set mesh_nodes($i,1) [expr $mesh_nodes($i,1)*$stretch_Y]
		set mesh_nodes($i,2) [expr $mesh_nodes($i,2)*$stretch_Z]
	}



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
	#for {set i 0} {$i < $mesh_nnodes} {incr i} {
		#set xx [discard_epsilon [expr $rotation(0,0)*$mesh_nodes($i,0) + $rotation(0,1)*$mesh_nodes($i,1) + $rotation(0,2)*$mesh_nodes($i,2)]]
		#set yy [discard_epsilon [expr $rotation(1,0)*$mesh_nodes($i,0) + $rotation(1,1)*$mesh_nodes($i,1) + $rotation(1,2)*$mesh_nodes($i,2)]]
		#set zz [discard_epsilon [expr $rotation(2,0)*$mesh_nodes($i,0) + $rotation(2,1)*$mesh_nodes($i,1) + $rotation(2,2)*$mesh_nodes($i,2)]]

		#set mesh_nodes($i,0) $xx;
		#set mesh_nodes($i,1) $yy;
		#set mesh_nodes($i,2) $zz;
	#}

	# rotation of nodes around X axis in YZ plane:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set yy [discard_epsilon [expr $ca*$mesh_nodes($i,1) - $sa*$mesh_nodes($i,2)]]
		set zz [discard_epsilon [expr $sa*$mesh_nodes($i,1) + $ca*$mesh_nodes($i,2)]]

		set mesh_nodes($i,1) $yy;
		set mesh_nodes($i,2) $zz;
	}

	# rotation of nodes around Y in ZX plane:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set zz [discard_epsilon [expr $cb*$mesh_nodes($i,2) - $sb*$mesh_nodes($i,0)]]
		set xx [discard_epsilon [expr $sb*$mesh_nodes($i,2) + $cb*$mesh_nodes($i,0)]]

		set mesh_nodes($i,2) $zz;
		set mesh_nodes($i,0) $xx;
	}

	# rotation of nodes around Z in XY plane:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set xx [discard_epsilon [expr $cc*$mesh_nodes($i,0) - $sc*$mesh_nodes($i,1)]]
		set yy [discard_epsilon [expr $sc*$mesh_nodes($i,0) + $cc*$mesh_nodes($i,1)]]

		set mesh_nodes($i,0) $xx;
		set mesh_nodes($i,1) $yy;
	}

	# setting the origin of the object-in-fluid:
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
		set mesh_nodes($i,0) [expr $mesh_nodes($i,0) + $origin_X]
		set mesh_nodes($i,1) [expr $mesh_nodes($i,1) + $origin_Y]
		set mesh_nodes($i,2) [expr $mesh_nodes($i,2) + $origin_Z]
	}

# Creating the list of edges	
	set mesh_nedges 0
	set mesh_nedges_bending 0
	
	for {set i 0} {$i < $mesh_ntriangles} {incr i} {
	    # Take a triangle and copy the nodes of the triangle to pa,pb,pc (point A, point B, point C)
	    set pa $mesh_triangles($i,0)
	    set pb $mesh_triangles($i,1)
	    set pc $mesh_triangles($i,2)
	    set is 0	

	    for {set j 0} {$j < $mesh_nedges} {incr j} {
			# Check if the edge AB or BA is in the current list of edges
			if {$mesh_edges($j,0) == $pa && $mesh_edges($j,1) == $pb} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pa
				set mesh_edges_bending($mesh_nedges_bending,1) $pb
				incr mesh_nedges_bending
			}
			if {$mesh_edges($j,1) == $pa && $mesh_edges($j,0) == $pb} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pa
				set mesh_edges_bending($mesh_nedges_bending,1) $pb
				incr mesh_nedges_bending
			}   
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
			if {$mesh_edges($j,0) == $pb && $mesh_edges($j,1) == $pc} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pb
				set mesh_edges_bending($mesh_nedges_bending,1) $pc
				incr mesh_nedges_bending
			}
			if {$mesh_edges($j,1) == $pb && $mesh_edges($j,0) == $pc} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pb
				set mesh_edges_bending($mesh_nedges_bending,1) $pc
				incr mesh_nedges_bending
			}  
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
			if {$mesh_edges($j,0) == $pa && $mesh_edges($j,1) == $pc} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pa
				set mesh_edges_bending($mesh_nedges_bending,1) $pc
				incr mesh_nedges_bending
			}
			if {$mesh_edges($j,1) == $pa && $mesh_edges($j,0) == $pc} { 
				set is 1
				set mesh_edges_bending($mesh_nedges_bending,0) $pa
				set mesh_edges_bending($mesh_nedges_bending,1) $pc
				incr mesh_nedges_bending
			}
		}

	    if {$is == 0} {     
		# If AC nor CA is in the list then add the edge AC to the list
		set mesh_edges($mesh_nedges,0) $pa
		set mesh_edges($mesh_nedges,1) $pc
		incr mesh_nedges
		}
	}    

puts "mesh_nedges = $mesh_nedges"
puts "mesh_nedges_bending = $mesh_nedges_bending"

	
#
#
#         ESPRESSO object creation
#	
# 	



# generating particles:
	if {$check_output == 1} { set f [open $createPart "w"] }
	set i $oif_firstPartId
	for {set i $oif_firstPartId} {$i < [expr $mesh_nnodes + $oif_firstPartId]} {incr i} {
		part $i pos [format %e [expr $mesh_nodes([expr $i - $oif_firstPartId],0)]] [format %e [expr $mesh_nodes([expr $i - $oif_firstPartId],1)]] [format %e [expr $mesh_nodes([expr $i - $oif_firstPartId],2)]] type $part_type mol $part_mol mass $part_mass
		if {$check_output == 1} { puts $f [format "part $i pos %e %e %e type $part_type mol $part_mol mass $part_mass" [expr $mesh_nodes([expr $i - $oif_firstPartId],0)] [expr $mesh_nodes([expr $i - $oif_firstPartId],1)] [expr $mesh_nodes([expr $i - $oif_firstPartId],2)]] }
		}  
	set oif_firstPartId [expr $oif_firstPartId + $mesh_nnodes]
	if {$check_output == 1} { close $f }

# generation of stretchlin force bonds
	# adapted by Cimo
	if { $kslin != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondSLIN "w"]
			set fpart [open $partSLIN "w"]
		}

		puts "generating stretchlin force bonds"
		set firstID_StrBond $oif_firstBondId
		# Stretchlin is coupled to the edges   
		set n_StrBond $mesh_nedges 
		set oif_firstBondId [expr $oif_firstBondId+$n_StrBond]
	
		set dist 0
		for {set i 0} {$i < $n_StrBond} {incr i} {
		    for {set k 0} {$k < 3} {incr k} {		
			set p1($k) [expr $mesh_nodes($mesh_edges($i,0),$k)]
			set p2($k) [expr $mesh_nodes($mesh_edges($i,1),$k)]
			}
		    # We need to compute the distance btw the vertices
		    set dist [expr sqrt(($p1(0)-$p2(0))*($p1(0)-$p2(0)) + ($p1(1)-$p2(1))*($p1(1)-$p2(1)) + ($p1(2)-$p2(2))*($p1(2)-$p2(2)))] 
		    set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
		    inter [expr $firstID_StrBond + $i] stretchlin_force [format %e $dist] [format %e $kslin]
		    #inter [expr $firstID_StrBond + $i] stretching_force [expr $dist] [expr $ks]
		    part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrBond + $i] [expr $mesh_edges($i,1) + $firstPartId]
		    if {$check_output == 1} { 
				puts $fbond [format "inter [expr $firstID_StrBond + $i] stretchlin_force %e %e" [expr $dist] [expr $kslin]]
				puts $fpart "part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrBond + $i] [expr $mesh_edges($i,1) + $firstPartId]"
			}
		}
		if {$check_output == 1} { 
			close $fpart
			close $fbond
		}
		
	}

		
# generation of stretching force bonds
	# implemented by Iveta
	if { $ks != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondS "w"]
			set fpart [open $partS "w"]
		}

		puts "generating stretching force bonds"
		set firstID_StrBond $oif_firstBondId
		# Stretching is coupled to the edges   
		set n_StrBond $mesh_nedges 
		set oif_firstBondId [expr $oif_firstBondId+$n_StrBond]
	
		set dist 0
		for {set i 0} {$i < $n_StrBond} {incr i} {
		    for {set k 0} {$k < 3} {incr k} {		
			set p1($k) [expr $mesh_nodes($mesh_edges($i,0),$k)]
			set p2($k) [expr $mesh_nodes($mesh_edges($i,1),$k)]
			}
		    # We need to compute the distance btw the vertices
		    set dist [expr sqrt(($p1(0)-$p2(0))*($p1(0)-$p2(0)) + ($p1(1)-$p2(1))*($p1(1)-$p2(1)) + ($p1(2)-$p2(2))*($p1(2)-$p2(2)))] 
		    set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
		    inter [expr $firstID_StrBond + $i] stretching_force [format %e $dist] [format %e $ks]
		    #inter [expr $firstID_StrBond + $i] stretching_force [expr $dist] [expr $ks]
		    part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrBond + $i] [expr $mesh_edges($i,1) + $firstPartId]
		    if {$check_output == 1} { 
				puts $fbond [format "inter [expr $firstID_StrBond + $i] stretching_force %e %e" [expr $dist] [expr $ks]]
				puts $fpart "part [expr $mesh_edges($i,0)+$firstPartId] bond [expr $firstID_StrBond + $i] [expr $mesh_edges($i,1) + $firstPartId]"
			}
		}
		if {$check_output == 1} { 
			close $fpart
			close $fbond
		}
		
	}

# generation of bending force bonds
	# implemented by Cimo
	if { $kb != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondB "w"]
			set fpart [open $partB "w"]
		}
		puts "generating bending force bonds"
		set firstID_BenBond $oif_firstBondId
		set n_BenBond $mesh_nedges_bending
			# Bending is coupled to the angles between triangles sharing the same edge
	    set oif_firstBondId [expr $oif_firstBondId + $n_BenBond]
	
		set phi 0.0
		for { set i 0} { $i < $n_BenBond} {incr i} { 
				#Run over all edges
			set p2id $mesh_edges_bending($i,0)
				#Put IDs of points to p2id,p3id
			set p3id $mesh_edges_bending($i,1) 
			for { set k 0} {$k < 3} {incr k} {	
				#Put coordinates of the edges's points
				set p2($k) $mesh_nodes($p2id,$k)
				set p3($k) $mesh_nodes($p3id,$k)
			}
			
			set detected 0
				#Number of detected triangles with current edge common
				# Algorithm is as follows: we run over all triangles and check whether two vertices are those from current edge. If we find such triangle, we put the ID of the third vertex to p1id and moreover we check if the orientation p1id, p2id p3id is the same as was in the triangle list (meaning, that we found one of the following three triples in the triangle list: p1id, p2id, p3id or p2id, p3id, p1id or p3id, p1id, p2id). If we have the same orientation, we set orient = 1, otherwise orient = -1.
				# Then we go further looking for the second triagle. The second triangle should have the opposite orientation.
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
	
			angle_btw_triangles P1 P2 P3 P4 phi

			# to replace lists with empty lists so they do not grow
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
			set P4 [lreplace $P4 0 2]

			inter [expr $firstID_BenBond + $i] bending_force [format %e $phi] [format %e $kb]
			set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
			part [expr $p2id + $firstPartId] bond [expr $firstID_BenBond + $i] [expr $p1id + $firstPartId] [expr $p3id + $firstPartId] [expr $p4id + $firstPartId]

			if {$check_output == 1} { 
				puts $fbond [format "inter [expr $firstID_BenBond + $i] bending_force %e %e" $phi $kb]
				puts $fpart "part [expr $p2id + $firstPartId] bond [expr $firstID_BenBond + $i] [expr $p1id + $firstPartId] [expr $p3id + $firstPartId] [expr $p4id + $firstPartId]"
			}
		}
		if {$check_output == 1} { 
			close $fpart
			close $fbond
		}
	}

# generation of local area force bonds
	# implemented by Cimo
	if {$kal != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondAlocal "w"]
			set fpart [open $partAlocal "w"]
		}
		puts "generating local area force bonds"
		set firstID_localAreaBond $oif_firstBondId
		set n_localAreaBond $mesh_ntriangles
			# Area is coupled to the triangles
	    set oif_firstBondId [expr $oif_firstBondId + $n_localAreaBond]
	
		set area 0.0
		for {set i 0} {$i < $n_localAreaBond} {incr i} {
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
			
			set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
			inter [expr $firstID_localAreaBond + $i] area_force_local $area $kal
			part [expr $mesh_triangles($i,0) + $firstPartId] bond [expr $firstID_localAreaBond + $i] [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]
			if {$check_output == 1} { 
				puts $fbond "inter [expr $firstID_localAreaBond + $i] area_force_local $area $kal"
				puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond [expr $firstID_localAreaBond + $i] [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]"
			}
		}
		if {$check_output == 1} {
			close $fbond
			close $fpart
		}
	}


# generation of global area force bonds
	# implemented by Cimo
	if {$kag != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondAglobal "w"]
			set fpart [open $partAglobal "w"]
		}
		puts "generating global area force bonds"
		set firstID_globalAreaBond $oif_firstBondId
		set n_globalAreaBond 1
	    set oif_firstBondId [expr $oif_firstBondId + $n_globalAreaBond]
	
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
		inter $firstID_globalAreaBond area_force_global $gl_area $kag
		if {$check_output == 1} { puts $fbond "inter $firstID_globalAreaBond area_force_global $gl_area $kag" }
			# First we need to set up the interaction and only afterward we can create the bonds
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
			if {$check_output == 1} { puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalAreaBond [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]" }
			part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalAreaBond [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]
		}
		if {$check_output == 1} { 
			close $fpart
			close $fbond
		}
	}
	
# generation of volume force bonds
	# implemented by Cimo
	if {$kv != 0.0} {
		if {$check_output == 1} { 
			set fbond [open $bondV "w"]
			set fpart [open $partV "w"]
		}
		puts "generating volume force bonds"
		set firstID_VolumeBond $oif_firstBondId
		set n_VolumeBond 1
		set oif_firstBondId [expr $oif_firstBondId + $n_VolumeBond]
	
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
		inter $firstID_VolumeBond volume_force $volume $kv
		if {$check_output == 1} { puts $fbond "inter $firstID_VolumeBond volume_force $volume $kv" }
			# First we need to set up the interaction and only afterward we can create the bonds
		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			set firstPartId [expr $oif_firstPartId - $mesh_nnodes]
			if {$check_output == 1} { puts $fpart "part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_VolumeBond [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) +$firstPartId]" }
			part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_VolumeBond [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) +$firstPartId]
		}
		if {$check_output == 1} { 
			close $fpart
			close $fbond
		}
	}
	
	
	# update global oif-variables
	incr oif_n_objects
	lappend oif_nnode $mesh_nnodes
	lappend oif_ntriangles $mesh_ntriangles
	lappend oif_nedges $mesh_nedges
}

