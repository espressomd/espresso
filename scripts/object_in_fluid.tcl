# Copyright (C) 2012,2013,2016 The ESPResSo project
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
    set oif_nparticles { } 
    # list of the numbers of nodes for each oif object

    global list oif_ntriangles
    set oif_ntriangles { }
    # list of the numbers of triangles for each oif object

    global list oif_nedges
    set oif_nedges { } 
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
    set oif_templates { }
    # 2D list of existing templates of oif objects containing rows with data. One row describes parameters of one object template. One row consists of X elements. The order of the elements is crucial. Structure of one row: num_of_particles, num_of_edges, num_of_triangles, ks, kslin, kb, kal, start_id_of_local_interaction, kag, kv, start_id_of_global_interaction, initial_surface, initial_volume, normal

    global list oif_objects
    set oif_objects { }
    # 2D list of existing objects containing rows with data. One row describes parameters of one object. One row consists of X elements. The order of the elements is crucial. Structure of one row: template_id, part_type, part_mass

    global list oif_template_nodes
    set oif_template_nodes { }
    # list containing triplets of coordinates of all nodes of all templates that will be used when creating objects. First num_of_nodes of them belong to the first template. Next num_of_nodes belong to the second template, etc. Note: these are not actual espresso particles, but rather node coordinates that will serve as template for the actual objects. 

    global list oif_template_edges
    set oif_template_edges { }
    # list containing pairs of particle IDs. The ordering works the same way as for particles. First num_of_edges edges belong to the first template, next num_of_edges belong to the second template, etc.

    global list oif_template_triangles
    set oif_template_triangles { }
    # list containing triplets of particle IDs. Again, in this list, triangles of all templates are arranged consequently and oif_template element num_of_triangles is used to determine which triangles belong to which template.

    global list oif_template_bending_incidences
    set oif_template_bending_incidences { }
    # a list of bending incidences (4 particle IDs in each row; number of rows equals number of rows in oif_template_edges). These are used so that during the creation of bending bonds one does not have to recalculate the bending incidences repeatedly. Number of bending interactions of a given template is equal to number of edges of this template.

    global list oif_template_3neighbors
    set oif_template_3neighbors { }
    # a list of 3 selected neighbors for each node (3 particle IDs in each row; number of rows equals number of rows in oif_template_nodes). These are used to compute the outward direction of object membrane. They are selected during creation of template, so that they do not have to be selected repeatedly for individual objects. Later they are used for cell-cell interactions.
    
    global list oif_object_starting_particles
    set oif_object_starting_particles { }
    # a list of particle indices denoting where each object's particles start in the whole list of particles

    global list oif_template_starting_triangles
    set oif_template_starting_triangles { }
    # a list of triangle indices denoting where each object's triangles start in the whole list of triangles

	global list oif_object_object_interactions
	set oif_object_object_interactions { }
	# a list indicating whether a non-bonded interaction has been set between objects. 
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
	global oif_template_nodes
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_template_3neighbors
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

        puts "oif_template_nodes"
	foreach item $oif_template_nodes {
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

	puts "oif_template_3neighbors"
	foreach item $oif_template_3neighbors {
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

proc norm {v} {
    # computes the norm of a 3-component vector
    upvar $v gv
    set vx [lindex $gv 0]
    set vy [lindex $gv 1]
    set vz [lindex $gv 2]
    set norm [expr sqrt($vx*$vx + $vy*$vy + $vz*$vz)]
    return $norm
}

proc distance {a b} {
    # computes the distance of two points
    upvar $a ga
    upvar $b gb
    
    set P1x [lindex $ga 0]
    set P1y [lindex $ga 1]
    set P1z [lindex $ga 2]
    set P2x [lindex $gb 0]
    set P2y [lindex $gb 1]
    set P2z [lindex $gb 2]

    set dist [expr sqrt(($P1x-$P2x)*($P1x-$P2x) + ($P1y-$P2y)*($P1y-$P2y) + ($P1z-$P2z)*($P1z-$P2z))]
    
    return $dist
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

#---------------------------------------------------------------------------
# calculation of stretching force

proc KS {lambda} {
    # Defined by (19) from Dupin2007
    set res [expr (pow($lambda,0.5) + pow($lambda,-2.5))/($lambda + pow($lambda,-3.))]
    return $res
}

proc calc_stretching_force { args } {
    # non-normalized stretching force
    set n_args 0
    foreach arg $args {
	incr n_args
    }
    if { $n_args != 9 } {
	puts "9 arguments are expected for calculation of stretching force"
	puts "ks, aX, aY, aZ, bX, bY, bZ, dist0, dist"
	return 0
    }

    set ks [lindex $args 0]
    set aX [lindex $args 1]
    set aY [lindex $args 2]
    set aZ [lindex $args 3]
    set bX [lindex $args 4]
    set bY [lindex $args 5]
    set bZ [lindex $args 6]
    set dist0 [lindex $args 7]
    set dist [lindex $args 8]
    
    set dr [expr $dist - $dist0]
    # two types of stretching:
    # nonlinear stretching:
    set lambda [expr 1.0*$dist/$dist0]
    set fac [expr (-$ks * [ KS $lambda ] * $dr )]
    # linear stretching: 
    #set fac [expr (-$ks * $dr)]
    
    set fX [expr ($fac*($bX-$aX)/$dist)]
    set fY [expr ($fac*($bY-$aY)/$dist)]
    set fZ [expr ($fac*($bZ-$aZ)/$dist)]
    
    set f [list $fX $fY $fZ]
    return $f
}
#---------------------------------------------------------------------------
# calculation of bending force

proc calc_bending_force { args } {
    set n_args 0
    foreach arg $args {
	incr n_args
    }
    if { $n_args != 15 } {
	puts "15 arguments are expected for calculation of bending force"
	puts "kb, aX, aY, aZ, bX, bY, bZ, cX, cY, cZ, dX, dY, dZ, phi0, phi"
	return 0
    }

    set kb [lindex $args 0]
    set aX [lindex $args 1]
    set aY [lindex $args 2]
    set aZ [lindex $args 3]
    set bX [lindex $args 4]
    set bY [lindex $args 5]
    set bZ [lindex $args 6]
    set cX [lindex $args 7]
    set cY [lindex $args 8]
    set cZ [lindex $args 9]
    set dX [lindex $args 10]
    set dY [lindex $args 11]
    set dZ [lindex $args 12]
    set phi0 [lindex $args 13]
    set phi [lindex $args 14]

    set pA [list $aX $aY $aZ]
    set pB [list $bX $bY $bZ]
    set pC [list $cX $cY $cZ]
    set pD [list $dX $dY $dZ]
    
    set n1 {0 0 0} 
    get_n_triangle pB pA pC n1
    set n1x [lindex $n1 0]
    set n1y [lindex $n1 1]
    set n1z [lindex $n1 2]
    set dn1 [norm n1]

    set n2 {0 0 0} 
    get_n_triangle pB pC pD n2
    set n2x [lindex $n2 0]
    set n2y [lindex $n2 1]
    set n2z [lindex $n2 2]
    set dn2 [norm n2]

    set angles [expr ($phi-$phi0)/$phi0]
    set fac [expr $kb * $angles]

    set f1x [expr $fac*$n1x/$dn1]
    set f1y [expr $fac*$n1y/$dn1]
    set f1z [expr $fac*$n1z/$dn1]
    set f2x [expr $fac*$n2x/$dn2]
    set f2y [expr $fac*$n2y/$dn2]
    set f2z [expr $fac*$n2z/$dn2]
    
    set f [list $f1x $f1y $f1z $f2x $f2y $f2z]
    return $f
}

#---------------------------------------------------------------------------
# calculation of local area force

proc calc_local_area_force { args } {
    set n_args 0
    foreach arg $args {
	incr n_args
    }
    if { $n_args != 12 } {
	puts "12 arguments are expected for calculation of local area force"
	puts "kal, aX, aY, aZ, bX, bY, bZ, cX, cY, cZ, A0, A"
	return 0
    }

    set kal [lindex $args 0]
    set a(0) [lindex $args 1]
    set a(1) [lindex $args 2]
    set a(2) [lindex $args 3]
    set b(0) [lindex $args 4]
    set b(1) [lindex $args 5]
    set b(2) [lindex $args 6]
    set c(0) [lindex $args 7]
    set c(1) [lindex $args 8]
    set c(2) [lindex $args 9]
    set A0 [lindex $args 10]
    set A [lindex $args 11]
    
    for {set i 0} {$i < 3} {incr i} {
	set h($i) [expr 1.0/3.0 *($a($i)+$b($i)+$c($i))]
    }
    
    set aa [expr ($A - $A0)/sqrt($A0)]

    # local area force for first node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$a($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f1($i) [expr $kal * $aa * $rh($i)/$hn]
    }

    # local area force for second node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$b($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f2($i) [expr $kal * $aa * $rh($i)/$hn]
    }

    # local area force for third node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$c($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f3($i) [expr $kal * $aa * $rh($i)/$hn]
    }
    
    set f [list $f1(0) $f1(1) $f1(2) $f2(0) $f2(1) $f2(2) $f3(0) $f3(1) $f3(2)]
    return $f
}

#---------------------------------------------------------------------------
# calculation of global area force

proc calc_global_area_force { args } {
    set n_args 0
    foreach arg $args {
	incr n_args
    }
    if { $n_args != 12 } {
	puts "12 arguments are expected for calculation of global area force"
	puts "kag, aX, aY, aZ, bX, bY, bZ, cX, cY, cZ, Ag0, Ag"
	return 0
    }

    set kag [lindex $args 0]
    set a(0) [lindex $args 1]
    set a(1) [lindex $args 2]
    set a(2) [lindex $args 3]
    set b(0) [lindex $args 4]
    set b(1) [lindex $args 5]
    set b(2) [lindex $args 6]
    set c(0) [lindex $args 7]
    set c(1) [lindex $args 8]
    set c(2) [lindex $args 9]
    set Ag0 [lindex $args 10]
    set Ag [lindex $args 11]
    
    for {set i 0} {$i < 3} {incr i} {
	set h($i) [expr 1.0/3.0 *($a($i)+$b($i)+$c($i))]
    }
    
    set aa [expr ($Ag - $Ag0)/$Ag0]

    # global area force for first node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$a($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f1($i) [expr $kag * $aa * $rh($i)/$hn]
    }

    # global area force for second node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$b($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f2($i) [expr $kag * $aa * $rh($i)/$hn]
    }

    # global area force for third node
    for {set i 0} {$i < 3} {incr i} {
	set rh($i) [expr $h($i)-$c($i)]
    }
    set rhv [list $rh(0) $rh(1) $rh(2)]
    set hn [norm rhv]
    for {set i 0} {$i < 3} {incr i} {
	set f3($i) [expr $kag * $aa * $rh($i)/$hn]
    }
    
    set f [list $f1(0) $f1(1) $f1(2) $f2(0) $f2(1) $f2(2) $f3(0) $f3(1) $f3(2)]
    return $f
}

#---------------------------------------------------------------------------
# calculation of volume force

proc calc_volume_force { args } {
    set n_args 0
    foreach arg $args {
	incr n_args
    }
    if { $n_args != 12 } {
	puts "12 arguments are expected for calculation of volume force"
	puts "kv, aX, aY, aZ, bX, bY, bZ, cX, cY, cZ, V0, V"
	return 0
    }

    set kv [lindex $args 0]
    set a(0) [lindex $args 1]
    set a(1) [lindex $args 2]
    set a(2) [lindex $args 3]
    set b(0) [lindex $args 4]
    set b(1) [lindex $args 5]
    set b(2) [lindex $args 6]
    set c(0) [lindex $args 7]
    set c(1) [lindex $args 8]
    set c(2) [lindex $args 9]
    set V0 [lindex $args 10]
    set V [lindex $args 11]
    
    set pA [list $a(0) $a(1) $a(2)]
    set pB [list $b(0) $b(1) $b(2)]
    set pC [list $c(0) $c(1) $c(2)]
    
    set n {0 0 0} 
    get_n_triangle pA pB pC n
    set nn(0) [lindex $n 0]
    set nn(1) [lindex $n 1]
    set nn(2) [lindex $n 2]
    set dn [norm n]

    set vv [expr ($V - $V0)/$V0]
    set A [area_triangle pA pB pC]    
   
    for {set i 0} {$i < 3} {incr i} {
	set f($i) [expr $kv * $vv * $A * $nn($i)/$dn / 3.0]
    }
    
    set fv [list $f(0) $f(1) $f(2)]
 
    return $fv
}

#-----------------------------------------------------------------------------
proc oif_create_template { args } {

# access global variables defined in oif_init
	global oif_first_inter_id
	global oif_n_templates
	global oif_templates
	global oif_template_nodes
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_template_3neighbors
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
	set ks_min 0.0
	set ks_max 0.0
	set ks_file ""
    set kslin 0.0
    set kslin_min 0.0
    set kslin_max 0.0
    set kslin_file ""
	set kb 0.0
	set kb_min 0.0
	set kb_max 0.0
	set kb_file ""
	set kal 0.0
	set kag 0.0
	set kv 0.0
	set check_output 0
	set template_id -1
	set normal -1

#--------------------------------------------------------------------------------------------
# reading the arguments. some of them are mandatory. we check for the mandatory arguments at the end of this section
	set pos 0
	while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"nodes-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing nodes-file"
					break
				}
				set filenamenodes [lindex $args $pos]
				incr pos
			}
			"ks" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value of stretching coefficient ks"
					break
				}
				set val [lindex $args $pos]
				if {[string is double $val] == 1} { 
					# uniform ks for all edges
					set ks $val
					incr pos
				} else { 
					set ks_file $val
					incr pos
					set ks_min [lindex $args $pos]
					incr pos
					set ks_max [lindex $args $pos]
					incr pos
					puts "ks_file $ks_file ks_min $ks_min ks_max $ks_max"
				}
			}
            "kslin" {
                incr pos
                if { $pos >= $n_args } {
                    puts "error in oif_create_template: missing value of linear stretching coefficient kslin"
                    break
                }
                set val [lindex $args $pos]
                if {[string is double $val] == 1} {
                    # uniform kslin for all edges
                    set kslin $val
                    incr pos
                } else {
                    set kslin_file $val
                    incr pos
                    set kslin_min [lindex $args $pos]
                    incr pos
                    set kslin_max [lindex $args $pos]
                    incr pos
                    puts "kslin_file $kslin_file kslin_min $kslin_min kslin_max $kslin_max"
                }
            }
			"kb" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value of bending coefficient kb"
					break
				}
				set val [lindex $args $pos]
				if {[string is double $val] == 1} { 
					# uniform kb for all edges
					set kb $val
					incr pos
				} else { 
					set kb_file $val
					incr pos
					set kb_min [lindex $args $pos]
					incr pos
					set kb_max [lindex $args $pos]
					incr pos
					puts "kb_file $kb_file kb_min $kb_min kb_max $kb_max"
				}
			}
			"kal" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value of local area coefficient kal"
					break
				}
				set kal [lindex $args $pos]
                incr pos
			}
			"kag" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value of global area coefficient kag"
					break
				}
				set kag [lindex $args $pos]
				incr pos
			}
			"kv" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value of volume coefficient kv"
					break
				}
				set kv [lindex $args $pos]
				incr pos
			}
			"triangles-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing triangles-file"
					break
				}
				set filenametriangles [lindex $args $pos]
				incr pos
			}
			"stretch" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: expecting 3 coefficients for stretching of the object"
					break
				}
				set stretch_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: expecting 3 coefficients for stretching of the object"
					break
				}
				set stretch_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: expecting 3 coefficients for stretching of the object"
					break
				}
				set stretch_Z [lindex $args $pos]
				incr pos
			}
		    "check" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing value for check_output"
					break
				}
				set check_output [lindex $args $pos]
				incr pos
			}
			"normal" {  
				incr pos
				set normal 1
			}
            "template-id" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_create_template: missing template id"
					break
				}
				set template_id [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_create_template: unknown keyword" 
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
		puts "error in oif_create_template: mandatory argument(s) missing"  
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
		set interLocal "TMP/interLocal[expr $template_id].txt"
		set interGlobal "TMP/interGlobal[expr $template_id].txt"
	}

#--------------------------------------------------------------------------------------------
# check: output all parameters:
	puts "The following oif_template will be created:"
	puts "	template-id: 	$template_id"
	puts "	nodes-file: 	$filenamenodes"
	puts "	triangles-file: $filenametriangles"
	puts "	ks:             $ks"
	puts "	kb:             $kb"
	puts "	kal:            $kal"
	puts "	kag:            $kag"
	puts "	kv:             $kv"
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
		        
	        list temp_node
	        lappend temp_node [expr $stretch_X*[lindex $line 0]]
		lappend temp_node [expr $stretch_Y*[lindex $line 1]]
		lappend temp_node [expr $stretch_Z*[lindex $line 2]]
	        
		# update global variable
		lappend oif_template_nodes $temp_node
	        
	        # clear list
		set temp_node [lreplace $temp_node 0 2]
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

#--------------------------------------------------------------------------------------------
	if { $normal != -1 } {

	# creating the list of neighbors
	# for each node we collect all its neighbors
	# and then we select three of them that are "best distributed around" the given node
	# we order them so that the normal of the triangle defined by them is pointing outward
	# these three are then saved in this order into the global list

	    puts "generating list of neighbors"
	    for {set i 0} {$i < $mesh_nnodes} {incr i} {
	    
	        set nneighbors 0
	        # cycle through edges and select those that contain node $i
	        for {set j 0} {$j < $mesh_nedges} {incr j} {
	        	# take an edge and copy the nodes of the edge to pa, pb
			set pa $mesh_edges($j,0)
			set pb $mesh_edges($j,1)
			# save neighbor if edge contains current node
			if {$pa == $i} {
			    set temp_neighbors($nneighbors) $pb
			    incr nneighbors
			}
			if {$pb == $i} {
			    set temp_neighbors($nneighbors) $pa
			    incr nneighbors
			}
	        }
	    
	        # create vectors to all neighbors and normalize them
	        for {set j 0} {$j < $nneighbors} {incr j} {
			for {set k 0} {$k < 3} {incr k} {
			    set vec_neighbor($j,$k) [expr $mesh_nodes($temp_neighbors($j),$k)-$mesh_nodes($i,$k)]
			}
			set temp_length [expr sqrt($vec_neighbor($j,0)*$vec_neighbor($j,0)+$vec_neighbor($j,1)*$vec_neighbor($j,1)+$vec_neighbor($j,2)*$vec_neighbor($j,2))]
			for {set k 0} {$k < 3} {incr k} {
			    set vec_neighbor($j,$k) [expr $vec_neighbor($j,$k)/$temp_length]
			}
	        }
	
	        # check all triplets of neighbors and select the one that is best spatially distributed
	        # by adding the corresponding three normalized vectors
	        # and selecting the one with smallest resultant vector
	        set best_neighbor(0) -1
	        set best_neighbor(1) -1
	        set best_neighbor(2) -1
	        set min_length 10000
	        for {set j 0} {$j < $nneighbors} {incr j} {
	        	for {set m [expr $j+1]} {$m < $nneighbors} {incr m} {
			    for {set n [expr $m+1]} {$n < $nneighbors} {incr n} {
				for {set k 0} {$k < 3} {incr k} {
				    set res($k) [expr $vec_neighbor($j,$k)+$vec_neighbor($m,$k)+$vec_neighbor($n,$k)]
				}
				set temp_length [expr sqrt($res(0)*$res(0)+$res(1)*$res(1)+$res(2)*$res(2))]
				if {$temp_length < $min_length} {
				    set min_length $temp_length
				    set best_neighbor(0) $temp_neighbors($j)
				    set best_neighbor(1) $temp_neighbors($m)
				    set best_neighbor(2) $temp_neighbors($n)
				}
			    }
			}
	        }
        
		# find one triangle that contains node $i and compute its normal vector - normal1
		for {set j 0} {$j < $mesh_ntriangles} {incr j} {
		    set pa $mesh_triangles($j,0)
		    set pb $mesh_triangles($j,1)
		    set pc $mesh_triangles($j,2)
		    if { $pa == $i || $pb == $i || $pc == $i} {
	                set normal1 [list 0.0 0.0 0.0]
		        set A [list $mesh_nodes($pa,0) $mesh_nodes($pa,1) $mesh_nodes($pa,2)]
		        set B [list $mesh_nodes($pb,0) $mesh_nodes($pb,1) $mesh_nodes($pb,2)]
		        set C [list $mesh_nodes($pc,0) $mesh_nodes($pc,1) $mesh_nodes($pc,2)]
		        get_n_triangle A B C normal1
		        set n1(0) [lindex $normal1 0]
		        set n1(1) [lindex $normal1 1]
		        set n1(2) [lindex $normal1 2]
		        set j $mesh_ntriangles
		    }
		}
	    
		# properly orient selected neighbors and update global list
		set normal2 [list 0.0 0.0 0.0]
		set A [list $mesh_nodes($best_neighbor(0),0) $mesh_nodes($best_neighbor(0),1) $mesh_nodes($best_neighbor(0),2)]
		set B [list $mesh_nodes($best_neighbor(1),0) $mesh_nodes($best_neighbor(1),1) $mesh_nodes($best_neighbor(1),2)]
		set C [list $mesh_nodes($best_neighbor(2),0) $mesh_nodes($best_neighbor(2),1) $mesh_nodes($best_neighbor(2),2)]
		get_n_triangle A B C normal2
		set n2(0) [lindex $normal2 0]
		set n2(1) [lindex $normal2 1]
		set n2(2) [lindex $normal2 2]
		set normn1 [expr sqrt($n1(0)*$n1(0) + $n1(1)*$n1(1) + $n1(2)*$n1(2))]
		set normn2 [expr sqrt($n2(0)*$n2(0) + $n2(1)*$n2(1) + $n2(2)*$n2(2))]
		set product [expr ($n1(0)*$n2(0)+$n1(1)*$n2(1)+$n1(2)*$n2(2))/($normn1 * $normn2)]
		set angle [expr acos($product)]
		list selected_neighbors
		lappend selected_neighbors $best_neighbor(0)
		set pi 3.1415926535897931
		if { $angle > $pi/2.0 } {
		    lappend selected_neighbors $best_neighbor(1)
		    lappend selected_neighbors $best_neighbor(2)
		} else {
		    lappend selected_neighbors $best_neighbor(2)
		    lappend selected_neighbors $best_neighbor(1)
		}
	        lappend oif_template_3neighbors $selected_neighbors
	    	        
	        # clean up
	        set selected_neighbors [lreplace $selected_neighbors 0 2]
	        for {set j 0} {$j < $nneighbors} {incr j} {
			set temp_neighbors($j) -1
	        }
	    }
	} else {
	    # note the nodes for which we do not want neighbors
	    list selected_neighbors
	    lappend selected_neighbors -1
	    lappend selected_neighbors -1
	    lappend selected_neighbors -1
	    for {set i 0} {$i < $mesh_nnodes} {incr i} {
		lappend oif_template_3neighbors $selected_neighbors
	    }
	    set selected_neighbors [lreplace $selected_neighbors 0 2]
	}

#----------------------------------------------------------------------------------------
# prepare template data for global variable oif_template 
	list template
	lappend template $mesh_nnodes
	lappend template $mesh_nedges
	lappend template $mesh_ntriangles
#-----------------------------------------------------------------------------------------
# generating oif_local_force interactions

    # only stretching and bending may be non-uniform
    if { $ks_file != "" } { lappend template 1.0 } else { lappend template $ks }
    if { $kslin_file != "" } { lappend template 1.0 } else { lappend template $kslin }
    if { $kb_file != "" } { lappend template 1.0 } else { lappend template $kb }
    lappend template $kal
    
    if { $ks != 0.0 && $kslin != 0.0} {
        puts "oif_create_template: warning: You have set both linear and non-linear stretching! Only one of the values ks and kslin should be non-zero."
    }
    
    if { $ks == 0.0 && $ks_file == "" && $kslin == 0.0 && $kslin_file == "" && $kb == 0.0 && $kb_file == "" && $kal == 0.0 } {
        set start_id_of_local_interaction -1
    } else {
        set start_id_of_local_interaction $oif_first_inter_id
    }
    lappend template $start_id_of_local_interaction
    
    if { $ks != 0.0 || $ks_file != "" || $kslin != 0.0 || $kslin_file != "" || $kb != 0.0 || $kb_file != "" || $kal != 0.0 } {
    
        if { $ks_file != ""} {
            # ks is nonuniform
            puts "reading nonuniform nonlinear stretching data"
            # reading file containing edges weights
            set fp [open $ks_file r]
            set file_data [read $fp]
            close $fp
            set data [split $file_data "\n"]
            set counter 0
            foreach line $data {
                if { [llength $line] == 1 } {
                    set ks_weight($counter) [lindex $line 0]
                    incr counter
                }
            }
            if { $counter != $mesh_nedges } {
                puts "Error in creating interactions for ks-nonuniform. Number of lines in $ks_file not equal to number of edges. Exiting."
                exit
            }
        }
    
        if { $kslin_file != ""} {
            # kslin is nonuniform
            puts "reading nonuniform linear stretching data"
            # reading file containing edges weights
            set fp [open $kslin_file r]
            set file_data [read $fp]
            close $fp
            set data [split $file_data "\n"]
            set counter 0
            foreach line $data {
                if { [llength $line] == 1 } {
                    set kslin_weight($counter) [lindex $line 0]
                    incr counter
                }
            }
            if { $counter != $mesh_nedges } {
                puts "Error in creating interactions for kslin-nonuniform. Number of lines in $kslin_file not equal to number of edges. Exiting."
                exit
            }
        }
    
        if { $kb_file != ""} {
            # kb is nonuniform
            puts "reading nonuniform bending data"
            # reading file containing edges weights
            set fp [open $kb_file r]
            set file_data [read $fp]
            close $fp
            set data [split $file_data "\n"]
            set counter 0
            foreach line $data {
                if { [llength $line] == 1 } {
                    set kb_weight($counter) [lindex $line 0]
                    incr counter
                }
            }
            if { $counter != $mesh_nedges } {
                puts "Error in creating interactions for kb-nonuniform. Number of lines in $kb_file not equal to number of edges. Exiting."
                exit
            }
        }
    
        puts "generating local force interactions"
        set n_LocalInter $mesh_nedges

        # Local interactions are coupled to the edges
        set oif_first_inter_id [expr $oif_first_inter_id + $n_LocalInter]
    
        if {$check_output == 1} {
            set out_file [open $interLocal "w"]
        }
    
        set phi 0.0
        for { set i 0} { $i < $n_LocalInter} {incr i} {
            # Run over all edges
            set p2id $mesh_edges($i,0)
            # Put IDs of points to p2id,p3id
            set p3id $mesh_edges($i,1)

            for { set k 0} {$k < 3} {incr k} {
                # Put coordinates of the edges' points
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
                        # if no triangle was detected
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
            
            # update global variable
            lappend oif_template_bending_incidences $temp_incidence
            # clear list
            set temp_incidence [lreplace $temp_incidence 0 3]
        
            angle_btw_triangles P1 P2 P3 P4 phi
            set area1 [area_triangle P1 P2 P3]
            set area2 [area_triangle P4 P2 P3]
            
            set dist [distance P2 P3]
                
            # to replace lists with empty lists so they do not grow
            set P1 [lreplace $P1 0 2]
            set P2 [lreplace $P2 0 2]
            set P3 [lreplace $P3 0 2]
            set P4 [lreplace $P4 0 2]
                
            if { $ks_file != ""} { set ks [expr $ks_min*(1 - $ks_weight($i)) + $ks_max*$ks_weight($i)] }
            if { $kslin_file != ""} { set kslin [expr $kslin_min*(1 - $kslin_weight($i)) + $kslin_max*$kslin_weight($i)] }
            if { $kb_file != ""} { set kb [expr $kb_min*(1 - $kb_weight($i)) + $kb_max*$kb_weight($i)] }

            inter [expr $start_id_of_local_interaction + $i] oif_local_force [format %e $dist] [format %e $ks] [format %e $kslin] [format %e $phi] [format %e $kb] [format %e $area1] [format %e $area2] [format %e $kal]
            if {$check_output == 1} {
                puts $out_file "inter [expr $start_id_of_local_interaction + $i] oif_local_force [format %e $dist] [format %e $ks] [format %e $kslin] [format %e $phi] [format %e $kb] [format %e $area1] [format %e $area2] [format %e $kal]"
            }
        }
        if {$check_output == 1} {
            close $out_file
        }
    }

#--------------------------------------------------------------------------------------------
# generation of oif_global_force interactions
	lappend template $kag
	lappend template $kv
	if { $kag == 0.0 && $kv == 0.0} {
	    set start_id_of_global_interaction -1
	    set gl_area -1
	    set gl_volume -1
	} else {
	    set start_id_of_global_interaction $oif_first_inter_id
	}
	lappend template $start_id_of_global_interaction

	if {$kag != 0.0 || $kv != 0.0} {
	    
	    puts "generating global force interactions"
	    set n_globalInter 1
	    set oif_first_inter_id [expr $oif_first_inter_id + $n_globalInter]

	    set area 0.0
	    set gl_area 0.0

	    set VOL_volume 0.0 
	    set VOL_hz 0.0
	    list VOL_norm
	    set VOL_dn 0.0
	    set VOL_drmax 0.0
	    set VOL_drtemp 0.0


		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
			for {set k 0} {$k < 3} {incr k} {
			    list P1
			    lappend P1 $mesh_nodes($mesh_triangles($i,0),$k)
			    list P2
			    lappend P2 $mesh_nodes($mesh_triangles($i,1),$k)
			    list P3
			    lappend P3 $mesh_nodes($mesh_triangles($i,2),$k)
			}
			
		# from volume
			set P10 [lindex $P1 0]
			set P11 [lindex $P1 1]
			set P12 [lindex $P1 2]
			set P20 [lindex $P2 0]
			set P21 [lindex $P2 1]
			set P22 [lindex $P2 2]
			set P30 [lindex $P3 0]
			set P31 [lindex $P3 1]
			set P32 [lindex $P3 2]
			set VOL_drtemp [expr sqrt(($P10 - $P20)*($P10 - $P20) + ($P11 - $P21)*($P11 - $P21) + ($P12 - $P22)*($P12 - $P22))]
			# distance P1P2
			if {$VOL_drmax < $VOL_drtemp} { set VOL_drmax $VOL_drtemp }
			set VOL_drtemp [expr sqrt(($P10 - $P30)*($P10 - $P30) + ($P11 - $P31)*($P11 - $P31) + ($P12 - $P32)*($P12 - $P32))]
			# distance P1P3
			if {$VOL_drmax < $VOL_drtemp} { set VOL_drmax $VOL_drtemp }
		#end from volume
			
			set area [area_triangle P1 P2 P3]
			set gl_area [expr $gl_area + $area]
		# from volume
			get_n_triangle P1 P2 P3 VOL_norm
			set VOL_norm0 [lindex $VOL_norm 0]
			set VOL_norm1 [lindex $VOL_norm 1]
			set VOL_norm2 [lindex $VOL_norm 2]
			set VOL_dn [expr sqrt($VOL_norm0*$VOL_norm0 + $VOL_norm1*$VOL_norm1 + $VOL_norm2*$VOL_norm2)]
			set VOL_hz [expr 1.0/3.0*($P12 + $P22 + $P32)]
			set VOL_volume [expr $VOL_volume + $area*$VOL_norm2/$VOL_dn*$VOL_hz]
            set gl_volume $VOL_volume
		# from volume
			# to replace lists with empty lists so they do not grow
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]
	    }

	    inter $start_id_of_global_interaction oif_global_forces $gl_area $kag $VOL_volume $kv

	    if {$check_output == 1} {
            set out_file [open $interGlobal "a"]
            puts $out_file "inter $start_id_of_global_interaction oif_global_forces $gl_area $kag $VOL_volume $kv"
            close $out_file
	    }
	    
	}

#--------------------------------------------------------------------------------------------
# update global oif-variables
	lappend template $gl_area
	lappend template $gl_volume
	if {$normal == -1} {
	    lappend template 0
	} else {
	    lappend template 1
	}
	lappend oif_templates $template
	incr oif_n_templates
} 

#--------------------------------------------------------------------------------------------
proc oif_add_object { args } {

# access global variables
	global oif_n_objects
	global oif_n_templates
	global oif_first_bond_id
	global oif_first_part_id
	global oif_first_inter_id
	global oif_templates
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_template_nodes
	global oif_template_edges
	global oif_template_triangles
	global oif_template_bending_incidences
	global oif_template_3neighbors
	global oif_object_starting_particles
	global oif_objects

#--------------------------------------------------------------------------------------------
# counts the number of arguments and initialize
	set n_args 0
	foreach arg $args {
		incr n_args
	}
	if { $n_args == 0 } {
		puts "Mandatory arguments are object-id, object template, origin, particle type"
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
					puts "error in oif_add_object: missing particle mass"
					break
				}
				set part_mass [lindex $args $pos]
				incr pos
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: missing object-id"
					break
				}
				set object_id [lindex $args $pos]
				set part_mol $object_id
				incr pos
			}
			"part-type" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: missing part-type"
					break
				}
				set part_type [lindex $args $pos]
				incr pos
			}
			"template-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: missing template-id"
					break
				}
				set template_id [lindex $args $pos]
				incr pos
			}
			"rotate" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: rotation input incorrect - 3 angles in radians expected"
					break
				}
				set rotate_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: rotation input incorrect - 3 angles in radians expected"
					break
				}
				set rotate_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: rotation input incorrect - 3 angles in radians expected"
					break
				}
				set rotate_Z [lindex $args $pos]
				incr pos
			}
			"origin" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: origin input incorrect - 3 coordinates expected"
					break
				}
				set origin_X [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: origin input incorrect - 3 coordinates expected"
					break
				}
				set origin_Y [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: origin input incorrect - 3 coordinates expected"
					break
				}
				set origin_Z [lindex $args $pos]
				incr pos
			}
			"check" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_add_object: output file name expected for check"
					break
				}
				set check_output [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_add_object: unknown keyword" 
				set pos $n_args
			}
		}  
	}


#--------------------------------------------------------------------------------------------
# checking wheter all mandatory arguments have been given
	set mandatory 1
	if { $origin_X == 0 && $origin_Y == 0 && $origin_Z == 0 } { set mandatory 0 }
	if { $part_type == "-1" } { set mandatory 0 }
	if { $template_id == "-1" } { set mandatory 0 }
	if { $object_id == "-1" } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_add_object: mandatory argument(s) missing" 
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
		set part "TMP/part$object_id.txt"
		set bondLocal "TMP/bondLocal$object_id.txt"
		set bondGlobal "TMP/bondGlobal$object_id.txt"
	}

#--------------------------------------------------------------------------------------------
# Check: output all parameters:
	puts "The following oif_object will be created:"
	puts "	object_id       $object_id"
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

 	# there are 14 elements in one row of oif_templates:
	# nnodes, nedges, ntriangles, ks, kslin, kb, kal, start_id_of_local_inter, kag, kv, start_id_of_global_inter, S_0, V_0, normal
	# get the data form oif_templates list
	set template [lindex $oif_templates $template_id]

	set mesh_nnodes [lindex $template 0]
	set mesh_nedges [lindex $template 1]
	set mesh_ntriangles [lindex $template 2]

	set start_id_of_local_interactions [lindex $template 7]
	set start_id_of_global_interactions [lindex $template 10]

	set start_id_of_particles 0
	set start_id_of_edges 0
	set start_id_of_triangles 0
	set start_id_of_bending_incidences 0
	for {set i 0} {$i < $template_id} {incr i} {
	    set start_id_of_particles [expr $start_id_of_particles + [lindex [lindex $oif_templates $i] 0]]
	    set start_id_of_edges [expr $start_id_of_edges + [lindex [lindex $oif_templates $i] 1]]
	    set start_id_of_triangles [expr $start_id_of_triangles + [lindex [lindex $oif_templates $i] 2]]
		set start_id_of_bending_incidences [expr $start_id_of_bending_incidences + [lindex [lindex $oif_templates $i] 1]]
	}
	 
 	# recover particles of the given template  
	for {set i 0} {$i < $mesh_nnodes} {incr i} {
	    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_particles+$i]]
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
    for {set i 0} {$i < $mesh_nedges} {incr i} {
        set bending_quartet [lindex $oif_template_bending_incidences [expr $start_id_of_bending_incidences+$i]]
		set bending_incidences($i,0) [lindex $bending_quartet 0]
		set bending_incidences($i,1) [lindex $bending_quartet 1]
		set bending_incidences($i,2) [lindex $bending_quartet 2]
		set bending_incidences($i,3) [lindex $bending_quartet 3]
	}
#--------------------------------------------------------------------------------------------
# determination of indices of mesh points that are extremal. Basically, this serves for approximate position of an object

	set xmin 100000
	set xmax -100000
	set ymin 100000
	set ymax -100000
	set zmin 100000
	set zmax -100000
	set xminID -1
	set xmaxID -1
	set yminID -1
	set ymaxID -1
	set zminID -1
	set zmaxID -1
	
	for {set ii 0 } { $ii < $mesh_nnodes } { incr ii } {
		if { $mesh_nodes($ii,0) > $xmax } { 
			set xmax $mesh_nodes($ii,0) 
			set xmaxID $ii
		}
		if { $mesh_nodes($ii,0) < $xmin } { 
			set xmin $mesh_nodes($ii,0) 
			set xminID $ii
		}
		if { $mesh_nodes($ii,1) > $ymax } { 
			set ymax $mesh_nodes($ii,1) 
			set ymaxID $ii
		}
		if { $mesh_nodes($ii,1) < $ymin } { 
			set ymin $mesh_nodes($ii,1) 
			set yminID $ii
		}
		if { $mesh_nodes($ii,2) > $zmax } { 
			set zmax $mesh_nodes($ii,2) 
			set zmaxID $ii
		}
		if { $mesh_nodes($ii,2) < $zmin } { 
			set zmin $mesh_nodes($ii,2) 
			set zminID $ii
		}
	}
	lappend object_data [list $xminID $xmaxID $yminID $ymaxID $zminID $zmaxID]

# recover list of neighbors of the given template
for {set i 0} {$i < $mesh_nnodes} {incr i} {
    set neighbors_triplet [lindex $oif_template_3neighbors [expr $start_id_of_particles+$i]]
    set mesh_3neighbors($i,0) [lindex $neighbors_triplet 0]
    set mesh_3neighbors($i,1) [lindex $neighbors_triplet 1]
    set mesh_3neighbors($i,2) [lindex $neighbors_triplet 2]
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
	set rotation(1,1) [expr $ca * $cc + $sa * $sb * $sc ]
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

#--------------------------------------------------------------------------------------------  
#
#
#         ESPRESSO object creation
#
# 
#--------------------------------------------------------------------------------------------
# generating particles:
	puts "generating particles"
	if {$check_output == 1} { set out_file [open $part "w"] }
	set i $oif_first_part_id

	# remember where this object's particles start:
	lappend oif_object_starting_particles $oif_first_part_id

	for {set i $oif_first_part_id} {$i < [expr $mesh_nnodes + $oif_first_part_id]} {incr i} {
		part $i pos [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],0)]] [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],1)]] [format %e [expr $mesh_nodes([expr $i - $oif_first_part_id],2)]] type $part_type mol $part_mol mass $part_mass
		if {$check_output == 1} { puts $out_file [format "part $i pos %e %e %e type $part_type mol $part_mol mass $part_mass" [expr $mesh_nodes([expr $i - $oif_first_part_id],0)] [expr $mesh_nodes([expr $i - $oif_first_part_id],1)] [expr $mesh_nodes([expr $i - $oif_first_part_id],2)]] }
		}  
	set oif_first_part_id [expr $oif_first_part_id + $mesh_nnodes]
	if {$check_output == 1} { close $out_file }

#--------------------------------------------------------------------------------------------
# generating local force bonds
	# template data
	# nnodes, nedges, ntriangles, ks, kslin, kb, kal, start_id_of_local_inter, kag, kv, start_id_of_global_inter, S_0, V_0, normal
	set ks_from_template [lindex $template 3]
    set kslin_from_template [lindex $template 4]
    set kb_from_template [lindex $template 5]
    set kal_from_template [lindex $template 6]
	if { $ks_from_template != 0.0 || $kslin_from_template != 0.0 || $kb_from_template != 0.0 || $kal_from_template != 0.0 } {
		
        if {$check_output == 1} {
		       	set out_file [open $bondLocal "w"]
		}
	   	
        puts "generating local force bonds"
		set firstID_LocalInter $start_id_of_local_interactions
        if { $firstID_LocalInter == -1 } { puts "oif_create_object: something went wrong with creating local interactions"; exit }
		# Local interactions are coupled to the edges
		set oif_first_bond_id [expr $oif_first_bond_id+$mesh_nedges]
		  
	   	for {set i 0} {$i < $mesh_nedges} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
            
            part [expr $bending_incidences($i,1) + $firstPartId] bond [expr $firstID_LocalInter + $i] [expr $bending_incidences($i,0) + $firstPartId] [expr $bending_incidences($i,2) + $firstPartId] [expr $bending_incidences($i,3) + $firstPartId]
		    if {$check_output == 1} { 
                puts $out_file "part [expr $bending_incidences($i,1) + $firstPartId] bond [expr $firstID_LocalInter + $i] [expr $bending_incidences($i,0) + $firstPartId] [expr $bending_incidences($i,2) + $firstPartId] [expr $bending_incidences($i,3) + $firstPartId]"
		    }
	   	}
		if {$check_output == 1} { 
		   close $out_file
	    }
	}

#--------------------------------------------------------------------------------------------
# generating global force bonds
	# template data
	# nnodes, nedges, ntriangles, ks, kslin, kb, kal, start_id_of_local_inter, kag, kv, start_id_of_global_inter, S_0, V_0, normal
	set kag_from_template [lindex $template 8]
	set kv_from_template [lindex $template 9]
	if { $kag_from_template != 0.0 || $kv_from_template != 0.0} {

		if {$check_output == 1} { 
		    set out_file [open $bondGlobal "w"]
	       	}
		
        puts "generating global force bonds"
		set firstID_globalInter $start_id_of_global_interactions
        if { $firstID_globalInter == -1 } { puts "oif_create_object: something went wrong with creating global interactions"; exit }
		set n_globalBond 1
		set oif_first_bond_id [expr $oif_first_bond_id + $n_globalBond]

		for {set i 0} {$i < $mesh_ntriangles} {incr i} {
		    set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
		    if {$check_output == 1} { 
			puts $out_file "part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]"
		    }
		    part [expr $mesh_triangles($i,0) + $firstPartId] bond $firstID_globalInter [expr $mesh_triangles($i,1) + $firstPartId] [expr $mesh_triangles($i,2) + $firstPartId]
		}
		if {$check_output == 1} { 
		    close $out_file
        }
	}
#--------------------------------------------------------------------------------------------
# computation of local outward direction of the membrane using the out_direction interaction and list of neighbors
    set normal_from_template [lindex $template 13]
    if { $normal_from_template == 1 } {
        puts "calculating outward directions"
        inter $oif_first_inter_id oif_out_direction
        
        set firstPartId [expr $oif_first_part_id - $mesh_nnodes]
        for {set i 0} {$i < $mesh_nnodes} {incr i} {
            part [expr $i + $firstPartId] bond $oif_first_inter_id [expr $mesh_3neighbors($i,0) + $firstPartId] [expr $mesh_3neighbors($i,1) + $firstPartId] [expr $mesh_3neighbors($i,2) + $firstPartId]
        }
        set oif_first_bond_id [expr $oif_first_bond_id + $mesh_nnodes]
        incr oif_first_inter_id
    }
#--------------------------------------------------------------------------------------------
# update global oif-variables
	incr oif_n_objects
	lappend oif_nparticles $mesh_nnodes
	lappend oif_ntriangles $mesh_ntriangles
	lappend oif_nedges $mesh_nedges
	lappend oif_objects $object_data
}

proc oif_template_change { args } {
	# access global variables
	global oif_templates
	
	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
	}
	if { $n_args == 0 } {
		puts "Mandatory argument is: template ID"
		return 0
	}

	set template_id -1
	set surface_perc 0
	set ks_perc 0
	set kal_perc 0
	set kag_perc 0
	set kv_perc 0
	
	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
	set pos 0
	while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"surface" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects a value for change of surface"
					break
				}
				set surface_perc [lindex $args $pos]
				incr pos
			}
			"ks-length" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects a value for change of stretching interactions"
					break
				}
				set ks_perc [lindex $args $pos]
				incr pos
			}
			"kal-area" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects a value for change of local area interactions"
					break
				}
				set kal_perc [lindex $args $pos]
				incr pos
			}
			"kag-area" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects a value for change of global area interaction"
					break
				}
				set kag_perc [lindex $args $pos]
				incr pos
			}
			"kv-volume" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects a value for change of volume interaction"
					break
				}
				set kv_perc [lindex $args $pos]
				incr pos
			}
			"template-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error: oif_template_change expects template-id"
					break
				}
				set template_id [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_template_change: unknown keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $template_id == -1 } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_template_change: mandatory argument(s) missing" 
		return
	}
	
	set template [lindex $oif_templates $template_id]
	
	# changing the stretching interactions for the given template
	if { $ks_perc != 0 } {
# template data
# nnodes, nedges, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter, S_0, V_0
	    # change cannot be too large (biological cells break at ~4% surface change) 
	    if { [expr abs($ks_perc) > 2] } {
		puts "oif_template_change: warning: change of stretching interactions might be too large"
	    }
	    set nedges [lindex $template 1]
	    set start_id_of_ks_interaction [lindex $template 4]
	    puts "	changing stretching interactions of template $template_id"
	    for {set i 0} {$i < $nedges} {incr i} {
		# read existing stretching interaction
		set old_inter [inter [expr $start_id_of_ks_interaction + $i]]
		set ks [lindex $old_inter 3]
		set old_dist [lindex $old_inter 2]
		# compute new relaxed length
		set dist [expr (1+0.01*$ks_perc)*$old_dist]
		# write new stretching interaction
		inter [expr $start_id_of_ks_interaction + $i] stretching_force [format %e $dist] [format %e $ks]
	    }
	}
	
	# changing the local area interactions for the given template
	if { $kal_perc != 0 } {
# template data
# nnodes, nedges, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter, S_0, V_0
	    if { [expr abs($kal_perc) > 4] } {
		puts "oif_template_change: warning: change of local area interactions might be too large"
	    }
	    set ntriangles [lindex $template 2]
	    set start_id_of_kal_interaction [lindex $template 8]
	    puts "	changing local area interactions of template $template_id"
	    for {set i 0} {$i < $ntriangles} {incr i} {
		# read existing stretching interaction
		set old_inter [inter [expr $start_id_of_kal_interaction + $i]]
		set kal [lindex $old_inter 3]
		set old_area [lindex $old_inter 2]
		# compute new relaxed area
		set area [expr (1+0.01*$kal_perc)*$old_area]
		# write new local area interaction
		inter [expr $start_id_of_kal_interaction + $i] area_force_local [format %e $area] [format %e $kal]
	    }
	}
	
	# changing the global area interaction for the given template
	if { $kag_perc != 0 } {
# template data
# nnodes, nedges, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter, S_0, V_0
	    if { [expr abs($kag_perc) > 4] } {
		puts "oif_template_change: warning: change of global area might be too large"
	    }
	    set kag [lindex $template 9]
	    set kag_inter_id [lindex $template 10]
	    set orig_surface [lindex $template 13]

	    puts "	changing global area interaction of template $template_id"
	    set new_surface [expr (1+0.01*$kag_perc)*$orig_surface]
	    inter $kag_inter_id area_force_global $new_surface $kag
	    lset oif_templates $template_id 13 $new_surface
	}
	
	# changing the volume for the given template
	if { $kv_perc != 0 } {
# template data
# nnodes, nedges, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter, S_0, V_0
	    if { [expr abs($kv_perc) > 8] } {
		puts "oif_template_change: warning: change of global volume might be too large"
	    }
	    set kv [lindex $template 11]
	    set kv_inter_id [lindex $template 12]
	    set orig_volume [lindex $template 14]
		
	    puts "	changing volume interaction of template $template_id"
	    set new_volume [expr (1+0.01*$kv_perc)*$orig_volume]
	    inter $kv_inter_id volume_force $new_volume $kv
	    lset oif_templates $template_id 14 $new_volume
	}
	
# changing the surface for the given template (includes change of stretching, local area and global area)
	if { $surface_perc != 0 } {
# template data
# nnodes, nedges, ntriangles, ks, start_id_of_ks_inter, kb, start_id_of_kb_inter, kal, start_id_of_kal_inter, kag, start_id_of_kag_inter, kv, start_id_of_kv_inter, S_0, V_0
	    set kag [lindex $template 9]
	    set kag_inter_id [lindex $template 10]
	    set orig_surface [lindex $template 13]
		
	    # surface change cannot be too large (biological cells break at ~4% surface change) 
	    if { [expr abs($surface_perc) > 3] } {
		puts "oif_template_change: warning: new surface should not differ by more than 3% from the old"
	    }
		
	    puts "	changing surface interaction of template $template_id"
	    set new_surface [expr (1+0.01*$surface_perc)*$orig_surface]
	    inter $kag_inter_id area_force_global $new_surface $kag
		
	    # also stretching and local area interactions will be changed
	    set stretching_factor [expr sqrt(1+0.01*$surface_perc)]
	    set nedges [lindex $template 1]
	    set start_id_of_ks_interaction [lindex $template 4]
	    puts "	changing stretching interactions of template $template_id"
	    for {set i 0} {$i < $nedges} {incr i} {
		# read existing stretching interaction
		set old_inter [inter [expr $start_id_of_ks_interaction + $i]]
		set ks [lindex $old_inter 3]
		set old_dist [lindex $old_inter 2]
		# compute new relaxed length
		set dist [expr $stretching_factor*$old_dist]
		# write new stretching interaction
		inter [expr $start_id_of_ks_interaction + $i] stretching_force [format %e $dist] [format %e $ks]
	    }
	    set ntriangles [lindex $template 2]
	    set start_id_of_kal_interaction [lindex $template 8]
	    puts "	changing local area interactions of template $template_id"
	    for {set i 0} {$i < $ntriangles} {incr i} {
		# read existing local area interaction
		set old_inter [inter [expr $start_id_of_kal_interaction + $i]]
		set kal [lindex $old_inter 3]
		set old_area [lindex $old_inter 2]
		# compute new relaxed local area
		set area [expr (1+0.01*$surface_perc)*$old_area]
		# write new local area interaction
		inter [expr $start_id_of_kal_interaction + $i] area_force_local [format %e $area] [format %e $kal]
	    } 
	}
}

proc oif_object_set { args } {
	# access global variables
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
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
		puts "Mandatory arguments are: object's ID"
		return 0
	}

	set objectID -1
	set force 0
	set velocity 0
	set mesh_nodes_file ""
	set origin 0
	set kill_motion 0
	set un_kill_motion 1
	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"force" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: force should be a 3d vector"
					break
				}
				set forceX [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: force should be a 3d vector"
					break
				}
				set forceY [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: force should be a 3d vector"
					break
				}
				set forceZ [lindex $args $pos]
				incr pos
				set force 1
			}
			"velocity" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: velocity should be a 3d vector"
					break
				}
				set velX [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: velocity should be a 3d vector"
					break
				}
				set velY [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: velocity should be a 3d vector"
					break
				}
				set velZ [lindex $args $pos]
				incr pos
				set velocity 1
			}
			"origin" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: origin should be a 3d vector"
					break
				}
				set originX [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: origin should be a 3d vector"
					break
				}
				set originY [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: origin should be a 3d vector"
					break
				}
				set originZ [lindex $args $pos]
				incr pos
				set origin 1
			}
			"mesh-nodes" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: mesh-nodes file missing"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"kill-motion" {  
				incr pos
				set kill_motion 1
			}
			"un-kill-motion" {  
				incr pos
				set un_kill_motion 1
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_set: object-id missing"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_object_set: incorrect keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_object_set: mandatory argument(s) missing" 
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

	if { $velocity == 1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			part $iii v $velX $velY $velZ
		}
	}

	if { $kill_motion == 1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			part $iii fix 1 1 1
		}
	}

	if { $un_kill_motion == 1} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			part $iii unfix
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
	# acces global variables:
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
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
	set vtk_pos_folded_file ""
	set vtk_aff_file ""
	set mesh_nodes_file ""

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"vtk-pos" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_output: missing output file for writing current positions"
					break
				}
				set vtk_pos_file [lindex $args $pos]
				incr pos
			}
			"vtk-pos-folded" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_output: missing output file for writing current positions"
					break
				}
				set vtk_pos_folded_file [lindex $args $pos]
				incr pos
			}
			"vtk-aff" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_output: missing output file for writing affinity bonds"
					break
				}
				set vtk_aff_file [lindex $args $pos]
				incr pos
			}
			"mesh-nodes" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_output: missing output file for writing mesh nodes"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_output: missing object-id"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_object_output: incorrect keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_object_output: missing mandatory arguments" 
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

	if { $vtk_pos_folded_file != ""} {
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
		set foldX [expr floor(1.0*[lindex $coords 0]/(1.0*[lindex [setmd box_l] 0]))]
		set foldY [expr floor(1.0*[lindex $coords 1]/(1.0*[lindex [setmd box_l] 1]))]
		set foldZ [expr floor(1.0*[lindex $coords 2]/(1.0*[lindex [setmd box_l] 2]))]
		# these gives how many times is the origin folded in all three directions

		set part [open $vtk_pos_folded_file "w"]
		puts $part "# vtk DataFile Version 3.0"
		puts $part "Data"
		puts $part "ASCII"
		puts $part "DATASET POLYDATA"
		puts $part "POINTS [lindex $oif_nparticles $objectID] float"
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set cpos [part $iii print pos]
			set posX [expr [lindex $cpos 0] - $foldX*[lindex [setmd box_l] 0]]
			set posY [expr [lindex $cpos 1] - $foldY*[lindex [setmd box_l] 1]]
			set posZ [expr [lindex $cpos 2] - $foldZ*[lindex [setmd box_l] 2]]
			puts $part "$posX $posY $posZ"
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


	if { $vtk_aff_file != ""} {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		set aff_nbonds 0

		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set tmppos [part $iii print pos]
			set tmpdata [part $iii print affinity]
			if { [lindex $tmpdata 0] != -1 } { 
				set aff_bondsA($aff_nbonds) $tmpdata 
				set aff_bondsB($aff_nbonds) $tmppos
				incr aff_nbonds
			}
		}
		
		if { $aff_nbonds > -1} {

			set part [open $vtk_aff_file "w"]
			puts $part "# vtk DataFile Version 3.0"
			puts $part "Data"
			puts $part "ASCII"
			puts $part "DATASET POLYDATA"
			puts $part "POINTS [expr 2*$aff_nbonds] float"
			for { set iii 0 } { $iii < $aff_nbonds } { incr iii } {
				puts $part "$aff_bondsA($iii)"
				puts $part "$aff_bondsB($iii)"
			}
			puts $part "LINES $aff_nbonds [expr 3*$aff_nbonds]"
			for { set iii 0 } { $iii < $aff_nbonds } { incr iii } {
				puts $part "2 [expr 2*$iii] [expr 2*$iii + 1]"
			}
			close $part
		}
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
	# access global variables 
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
	global oif_templates
	global oif_template_starting_triangles
	global oif_template_bending_incidences
	global oif_template_nodes
	global oif_template_triangles
	global oif_template_edges
	global timestep

	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "Mandatory arguments are: object's ID"
		return 0
	}

	set objectID -1
	set center_location 0
	set pos_bounds ""
	set volume -1
	set surface -1
	set elastic_forces ""
	set output_file ""
	set velocity 0
	set affinity ""
	set approx_pos -1
	set get_first_part_id 0
	set n_nodes 0
	set aff_type -1
	set aff_kappa -1
	set aff_r0 -1
	set aff_Kon -1
	set aff_K0 -1
	set aff_Fd -1
	set aff_file "aff-default.dat"

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"origin" {  
				incr pos
				set center_location 1
			}
			"approx-pos" {  
				incr pos
				set approx_pos 1
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
			"first-particle-id" {  
				incr pos
				set get_first_part_id 1
			}
			"n-nodes" {  
				incr pos
				set n_nodes 1
			}
			"affinity" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: type of output for affinity missing (possible choices: nbonds, aver-bond-length, all)"
					break
				}
				set affinity [lindex $args $pos]
				incr pos
			}
			"active-aff-bonds" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: relaxed length is required"
					break
				}
				set aff_type [lindex $args $pos]
				incr pos
				set aff_kappa [lindex $args $pos]
				incr pos
				set aff_r0 [lindex $args $pos]
				incr pos
				set aff_Kon [lindex $args $pos]
				incr pos
				set aff_K0 [lindex $args $pos]
				incr pos
				set aff_Fd [lindex $args $pos]
				incr pos
				set aff_file [lindex $args $pos]
				incr pos
			}
			"pos-bounds" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: type of pos-bounds is missing"
					puts "pos-bounds options: all, x-min, x-max, y-min, y-max, z-min, z-max"
					break
				}
				set pos_bounds [lindex $args $pos]
				incr pos
			}
			"elastic-forces" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: type of elastic-forces or f-metric is missing"
					puts "elastic-forces options: ks, kb, kal, kag, kv, total, ks-fmetric, kb-fmetric, kal-fmetric, kag-fmetric, kv-fmetric"
					break
				}
				set elastic_forces [lindex $args $pos]
				incr pos
				if { $elastic_forces == "ks" || $elastic_forces == "kb" || $elastic_forces == "kal" || \
				    $elastic_forces == "kag" || $elastic_forces == "kv" || $elastic_forces == "total"} {
				    if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: output file for elastic-forces is missing"
					break
				    }
				    set output_file [lindex $args $pos]
				    incr pos
				}
			}
			"object-id" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_object_analyze: object-id is missing"
					break
				}
				set objectID [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_object_analyze: incorrect keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $objectID == -1 } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_object_analyze: mandatory argument(s) missing" 
		return
	}

	if { $get_first_part_id == 1 } {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		return $firstPartId
	}

	if { $n_nodes == 1 } {
		set nnode [lindex $oif_nparticles $objectID]
		return $nnode
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

	if { $aff_type % 10 == 2 } {
		set fout [open $aff_file a]
		puts $fout "$aff_type $aff_kappa $aff_r0 $aff_Kon $aff_K0 $aff_Fd"
		set n_active_bonds 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set bond_site [part $iii print affinity]
			set bsx [lindex $bond_site 0]
			set bsy [lindex $bond_site 1]
			set bsz [lindex $bond_site 2]
			if { $bsx != -1 || $bsy != -1 || $bsz != -1} {

				set coords [part $iii print pos]
				set cx [lindex $coords 0]
				set cy [lindex $coords 1]
				set cz [lindex $coords 2]
				set bond_length [expr sqrt(($cx-$bsx)*($cx-$bsx) + ($cy-$bsy)*($cy-$bsy) + ($cz-$bsz)*($cz-$bsz))]
				if { $bond_length > $aff_r0 } {
					set tmpF [expr $aff_kappa*($bond_length - $aff_r0)]
					set tmpKoff [expr $aff_K0*exp( $tmpF / $aff_Fd)]
					set tmpPoff [expr 1.0 - exp( - $tmpKoff*$timestep)]
					puts $fout  "$iii: len $bond_length Koff $tmpKoff Poff $tmpPoff "
					incr n_active_bonds
				}
			}
		}
		close $fout
		return $n_active_bonds
	}

	if { $aff_type % 10 == 5 } {
		set fout [open $aff_file a]
		puts $fout "$aff_type $aff_kappa $aff_r0 $aff_Kon $aff_K0 $aff_Fd"
		set n_active_bonds 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set bond_site [part $iii print affinity]
			set bsx [lindex $bond_site 0]
			set bsy [lindex $bond_site 1]
			set bsz [lindex $bond_site 2]
			if { $bsx != -1 || $bsy != -1 || $bsz != -1} {

				set coords [part $iii print pos]
				set cx [lindex $coords 0]
				set cy [lindex $coords 1]
				set cz [lindex $coords 2]
				set bond_length [expr sqrt(($cx-$bsx)*($cx-$bsx) + ($cy-$bsy)*($cy-$bsy) + ($cz-$bsz)*($cz-$bsz))]
				if { $bond_length > [expr 0.75*$aff_r0] } {
					set tmpF [expr $aff_kappa*($bond_length - 0.75*$aff_r0)]
					set tmpKoff [expr $aff_K0*exp( $tmpF / $aff_Fd)]
					set tmpPoff [expr 1.0 - exp( - $tmpKoff*$timestep)]
					puts $fout "$iii: len $bond_length Koff $tmpKoff Poff $tmpPoff "
					incr n_active_bonds
				}
			}
		}
		close $fout
		return $n_active_bonds
	}

	if { $approx_pos == 1 } {
		set appX 0
		set appY 0
		set appZ 0
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set object_data [lindex $oif_objects $objectID] 
		set extremal_indices [lindex $object_data 3]
		set xminID [expr $firstPartId + [lindex $extremal_indices 0]]
		set xmaxID [expr $firstPartId + [lindex $extremal_indices 1]]
		set yminID [expr $firstPartId + [lindex $extremal_indices 2]]
		set ymaxID [expr $firstPartId + [lindex $extremal_indices 3]]
		set zminID [expr $firstPartId + [lindex $extremal_indices 4]]
		set zmaxID [expr $firstPartId + [lindex $extremal_indices 5]]

		set coords [part $xminID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set coords [part $xmaxID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set coords [part $yminID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set coords [part $ymaxID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set coords [part $zminID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set coords [part $zmaxID print pos]
		set appX [expr $appX + [lindex $coords 0]]
		set appY [expr $appY + [lindex $coords 1]]
		set appZ [expr $appZ + [lindex $coords 2]]

		set appX [expr $appX/6]
		set appY [expr $appY/6]
		set appZ [expr $appZ/6]
		set coords [list $appX $appY $appZ]
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

	if { $affinity != "" } {
		set firstPartId [lindex $oif_object_starting_particles $objectID]
		set nnode [lindex $oif_nparticles $objectID]
		set n_bonds 0
		set total_length 0
		set min_bond_length 10000000
		set max_bond_length -10000000
		set min_x_bond 10000000
		set max_x_bond -10000000
		for { set iii $firstPartId } { $iii < [expr $firstPartId + $nnode] } { incr iii } {
			set bond_site [part $iii print affinity]
			set bsiteX [lindex $bond_site 0]
			set bsiteY [lindex $bond_site 1]
			set bsiteZ [lindex $bond_site 2]
			if { $bsiteX != -1 || $bsiteY != -1 || $bsiteZ != -1 } {
				incr n_bonds
				set part_pos [part $iii print pos]
				set posX [lindex $part_pos 0]
				set posY [lindex $part_pos 1]
				set posZ [lindex $part_pos 2]
				set curr_bond_length [expr sqrt(($posX - $bsiteX)*($posX - $bsiteX) + ($posY - $bsiteY)*($posY - $bsiteY) + ($posZ - $bsiteZ)*($posZ - $bsiteZ))]
				set total_length [expr $total_length + $curr_bond_length] 
				if { $curr_bond_length < $min_bond_length } { set min_bond_length $curr_bond_length}
				if { $curr_bond_length > $max_bond_length } { set max_bond_length $curr_bond_length}
				if { $bsiteX < $min_x_bond } { set min_x_bond $bsiteX}
				if { $bsiteX > $max_x_bond } { set max_x_bond $bsiteX}
			}
		}
		if { $n_bonds > 0 } { set total_length [expr 1.0*$total_length / (1.0*$n_bonds)] }
		if { $affinity == "nbonds" } { return $n_bonds }
		if { $affinity == "aver-bond-length" } { 
			return $total_length 
		}
		if { $affinity == "min-bond-length" } { 
			return $min_bond_length 
		}
		if { $affinity == "max-bond-length" } { 
			return $max_bond_length 
		}
		if { $affinity == "min-x-bond" } { 
			return $min_x_bond 
		}
		if { $affinity == "max-x-bond" } { 
			return $max_x_bond 
		}
		if { $affinity == "all" } { 
			set answer [list "nbonds" $n_bonds "bond-length" $total_length "min-bond-length" $min_bond_length "max-bond-length" $max_bond_length "min-x-bond" $min_x_bond "max-x-bond" $max_x_bond]
			return $answer 
		}
		
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

	if { $elastic_forces != ""} {
		set result -1
		if { $elastic_forces == "ks" || $elastic_forces == "total" || $elastic_forces == "ks-fmetric" } {
			set result 0

			# recover data of this object
			set object [lindex $oif_objects $objectID]
			set template_id [lindex $object 0]
			set template [lindex $oif_templates $template_id]
			set nnodes [lindex $template 0]
			set nedges [lindex $template 1]
			set ks [lindex $template 3]

			set start_id_of_nodes 0
			set start_id_of_edges 0
			for {set i 0} {$i < $template_id} {incr i} {
			    set start_id_of_nodes [expr $start_id_of_nodes + [lindex [lindex $oif_templates $i] 0]]
			    set start_id_of_edges [expr $start_id_of_edges + [lindex [lindex $oif_templates $i] 1]]
			}

			# recover particles from this object's template  
			for {set i 0} {$i < $nnodes} {incr i} {
			    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_nodes+$i]]
			    set mesh_nodes($i,0) [lindex $node_triplet 0] 
			    set mesh_nodes($i,1) [lindex $node_triplet 1]
			    set mesh_nodes($i,2) [lindex $node_triplet 2]
			}

			# recover edges from this object's template
			for {set i 0} {$i < $nedges} {incr i} {
			    set edge_pair [lindex $oif_template_edges [expr $start_id_of_edges+$i]]
			    set mesh_edges($i,0) [lindex $edge_pair 0]
			    set mesh_edges($i,1) [lindex $edge_pair 1]
			}

			# initialize the mesh_stretching list
			for { set i 0 } { $i < $nnodes } { incr i } {
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_stretching($i,$j) 0.0
			    }
			}
			  
			set start_id_of_particles [lindex $oif_object_starting_particles $objectID]  
			    
			for { set i 0 } { $i < $nedges } { incr i } {
			    # Take an edge and copy the nodes of the edge to pA, pB (point A, point B)
			    set pA $mesh_edges($i,0)
			    set pB $mesh_edges($i,1)
			    
			    # Get the current position of the same two points
			    set currA [part [expr $start_id_of_particles + $pA] print pos]
			    set currAx [lindex $currA 0]
			    set currAy [lindex $currA 1]
			    set currAz [lindex $currA 2]
			    set currB [part [expr $start_id_of_particles + $pB] print pos]
			    set currBx [lindex $currB 0]
			    set currBy [lindex $currB 1]
			    set currBz [lindex $currB 2]

			    # Calculate stretching force between point A and point B
			    set orig_dist 0.0
			    for { set j 0 } { $j < 3 } { incr j } {
				set orig_dist [expr ($orig_dist + ($mesh_nodes($pB,$j) - $mesh_nodes($pA,$j))*($mesh_nodes($pB,$j) - $mesh_nodes($pA,$j)))]
			    }
			    set orig_dist [expr sqrt($orig_dist)]
			    
			    set curr_dist 0.0
			    set curr_dist [expr ($curr_dist + ($currBx - $currAx)*($currBx - $currAx))]
			    set curr_dist [expr ($curr_dist + ($currBy - $currAy)*($currBy - $currAy))]
			    set curr_dist [expr ($curr_dist + ($currBz - $currAz)*($currBz - $currAz))]
			    set curr_dist [expr sqrt($curr_dist)]
			    
			    set fTemp [calc_stretching_force $ks $currAx $currAy $currAz $currBx $currBy $currBz $orig_dist $curr_dist]
			    set f(0) [lindex $fTemp 0]
			    set f(1) [lindex $fTemp 1]
			    set f(2) [lindex $fTemp 2]
			      
			    # Save this stretching force into the array of stretching forces 
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_stretching($pA,$j) [expr $mesh_stretching($pA,$j) + $f($j)]
				set mesh_stretching($pB,$j) [expr $mesh_stretching($pB,$j) - $f($j)]
			    }
			    
			}

			# output here only if elastic forces are set to ks
			if { $elastic_forces == "ks" } {
			    # determine whether .vtk or other file type will be written
			    if { [ string first ".vtk" $output_file ] > -1 } {
			        # write stretching forces into a .vtk file for visualization
			        set out_f [open $output_file "w"]
			        puts $out_f "# vtk DataFile Version 3.0"
			        puts $out_f "Data"
			        puts $out_f "ASCII"
			        puts $out_f "DATASET POLYDATA"
			        puts $out_f "POINTS $nnodes float"
			        # writing points
			        for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
					puts $out_f "[part $i print pos]"
			        }
			        set ntriangles [lindex $template 2]
			        puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			        # extracts the starting position of the triangles for the template_id
			        set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			        # writing triangles
			        for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
					puts $out_f "3 [lindex $oif_template_triangles $i]"
			        }
			        # writing stretching forces
			        puts $out_f "POINT_DATA $nnodes"
			        puts $out_f "SCALARS sample_scalars float 1"
			        puts $out_f "LOOKUP_TABLE default"
			        for { set i 0 } { $i < $nnodes } { incr i } {
			        	set total_stretching 0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_stretching [expr $total_stretching + $mesh_stretching($i,$j)*$mesh_stretching($i,$j)]
					}
					set total_stretching [expr sqrt($total_stretching)]
					puts $out_f "$total_stretching"
				}
				close $out_f
			    } else {
			        # write stretching forces into a text file
			        set out_f [open $output_file "w"]
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_stretching 0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_stretching [expr $total_stretching + $mesh_stretching($i,$j)*$mesh_stretching($i,$j)]
					}
					set total_stretching [expr sqrt($total_stretching)]
					puts $out_f "$total_stretching"
				}
				close $out_f
			    }
			}
			if { $elastic_forces == "ks-fmetric" } {
			    set ks_fmetric 0.0 
			    for { set i 0 } { $i < $nnodes } { incr i } {
				set total_stretching 0.0
				for { set j 0 } { $j < 3 } { incr j } {
				    # compute the magnitude of the resultant force in one node 
				    set total_stretching [expr $total_stretching + $mesh_stretching($i,$j)*$mesh_stretching($i,$j)]
				}
				set total_stretching [expr sqrt($total_stretching)]
				set ks_fmetric [expr $ks_fmetric + $total_stretching]
			    }
			    set result $ks_fmetric
			}
		}
		if { $elastic_forces == "kb" || $elastic_forces == "total" || $elastic_forces == "kb-fmetric" } {
			set result 0

			# recover data of this object
			set object [lindex $oif_objects $objectID]
			set template_id [lindex $object 0]
			set template [lindex $oif_templates $template_id]
			set nnodes [lindex $template 0]
			set nedges [lindex $template 1]
			set kb [lindex $template 5]

			set start_id_of_nodes 0
			set start_id_of_edges 0
			set start_id_of_bending_incidences 0
			for {set i 0} {$i < $template_id} {incr i} {
			    set start_id_of_nodes [expr $start_id_of_nodes + [lindex [lindex $oif_templates $i] 0]]
			    set start_id_of_edges [expr $start_id_of_edges + [lindex [lindex $oif_templates $i] 1]]
			    if { [lindex [lindex $oif_templates $i] 5] != 0.0 } {
				set start_id_of_bending_incidences [expr $start_id_of_bending_incidences + [lindex [lindex $oif_templates $i] 1]]
			    }
			}

			# recover particles from this object's template  
			for {set i 0} {$i < $nnodes} {incr i} {
			    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_nodes+$i]]
			    set mesh_nodes($i,0) [lindex $node_triplet 0] 
			    set mesh_nodes($i,1) [lindex $node_triplet 1]
			    set mesh_nodes($i,2) [lindex $node_triplet 2]
			}

			# recover bending incidences from this object's template
			for {set i 0} {$i < $nedges} {incr i} {
			    set incidence [lindex $oif_template_bending_incidences [expr $start_id_of_bending_incidences + $i]]
			    set mesh_incidences($i,0) [lindex $incidence 0]
			    set mesh_incidences($i,1) [lindex $incidence 1]
			    set mesh_incidences($i,2) [lindex $incidence 2]
			    set mesh_incidences($i,3) [lindex $incidence 3]
			}

			# initialize the mesh_bending list
			for { set i 0 } { $i < $nnodes } { incr i } {
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_bending($i,$j) 0.0
			    }
			}
			  
			set start_id_of_particles [lindex $oif_object_starting_particles $objectID]    
			  
			set orig_phi 0.0
			set curr_phi 0.0
			for { set i 0 } { $i < $nedges } { incr i } {
			    # Take an edge and copy the bending incidences to pA, pB, pC, pD
			    set pA $mesh_incidences($i,0)
			    set pB $mesh_incidences($i,1)
			    set pC $mesh_incidences($i,2)
			    set pD $mesh_incidences($i,3)
			   
			    set origA [list $mesh_nodes($pA,0) $mesh_nodes($pA,1) $mesh_nodes($pA,2)]
			    set origB [list $mesh_nodes($pB,0) $mesh_nodes($pB,1) $mesh_nodes($pB,2)]
			    set origC [list $mesh_nodes($pC,0) $mesh_nodes($pC,1) $mesh_nodes($pC,2)]
			    set origD [list $mesh_nodes($pD,0) $mesh_nodes($pD,1) $mesh_nodes($pD,2)]
			    
			    # Get the current position of the same four points
			    set currA [part [expr $start_id_of_particles + $pA] print pos]
			    set currAx [lindex $currA 0]
			    set currAy [lindex $currA 1]
			    set currAz [lindex $currA 2]
			    set currB [part [expr $start_id_of_particles + $pB] print pos]
			    set currBx [lindex $currB 0]
			    set currBy [lindex $currB 1]
			    set currBz [lindex $currB 2]
			    set currC [part [expr $start_id_of_particles + $pC] print pos]
			    set currCx [lindex $currC 0]
			    set currCy [lindex $currC 1]
			    set currCz [lindex $currC 2]
			    set currD [part [expr $start_id_of_particles + $pD] print pos]
			    set currDx [lindex $currD 0]
			    set currDy [lindex $currD 1]
			    set currDz [lindex $currD 2]
			    
			    angle_btw_triangles origA origB origC origD orig_phi
			    angle_btw_triangles currA currB currC currD curr_phi
			    
			    # Calculate bending force among these four points
			    set fTemp [calc_bending_force $kb $currAx $currAy $currAz $currBx $currBy $currBz \
				       $currCx $currCy $currCz $currDx $currDy $currDz $orig_phi $curr_phi]
			    # first three components of fTemp are force on triangle 1, last three components are force on triangle 2
			    set f(0) [lindex $fTemp 0]
			    set f(1) [lindex $fTemp 1]
			    set f(2) [lindex $fTemp 2]
			    set f(3) [lindex $fTemp 3]
			    set f(4) [lindex $fTemp 4]
			    set f(5) [lindex $fTemp 5]
			    
			    # Save these bending forces into the array of bending forces 
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_bending($pA,$j) [expr $mesh_bending($pA,$j) + $f($j)]
				set mesh_bending($pB,$j) [expr $mesh_bending($pB,$j) + 0.5*$f($j) + 0.5*$f([expr $j+3])]
				set mesh_bending($pC,$j) [expr $mesh_bending($pC,$j) + 0.5*$f($j) + 0.5*$f([expr $j+3])]
				set mesh_bending($pD,$j) [expr $mesh_bending($pD,$j) + $f([expr $j+3])]
			    }
			    
			}

			# output here only if elastic forces are set to kb
			if { $elastic_forces == "kb" } {
			    # determine whether .vtk or other file type will be written
			    if { [ string first ".vtk" $output_file ] > -1 } {
			        # write bending forces into a .vtk file for visualization
			        set out_f [open $output_file "w"]
			        puts $out_f "# vtk DataFile Version 3.0"
			        puts $out_f "Data"
			        puts $out_f "ASCII"
			        puts $out_f "DATASET POLYDATA"
			        puts $out_f "POINTS $nnodes float"
			        # writing points
			        for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
					puts $out_f "[part $i print pos]"
			        }
			        set ntriangles [lindex $template 2]
			        puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			        # extracts the starting position of the triangles for the template_id
			        set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			        # writing triangles
			        for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
					puts $out_f "3 [lindex $oif_template_triangles $i]"
			        }
			        # writing bending forces
			        puts $out_f "POINT_DATA $nnodes"
			        puts $out_f "SCALARS sample_scalars float 1"
			        puts $out_f "LOOKUP_TABLE default"
			        for { set i 0 } { $i < $nnodes } { incr i } {
			        	set total_bending 0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_bending [expr $total_bending + $mesh_bending($i,$j)*$mesh_bending($i,$j)]
					}
					set total_bending [expr sqrt($total_bending)]
					puts $out_f "$total_bending"
			        }
			        close $out_f
			    } else {
			        # write bending forces into a text file
			        set out_f [open $output_file "w"]
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_bending 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_bending [expr $total_bending + $mesh_bending($i,$j)*$mesh_bending($i,$j)]
					}
					set total_bending [expr sqrt($total_bending)]
					puts $out_f "$total_bending"
			        }
			        close $out_f
			    }
			}
			if { $elastic_forces == "kb-fmetric" } {
			    set kb_fmetric 0.0 
			    for { set i 0 } { $i < $nnodes } { incr i } {
				set total_bending 0.0
				for { set j 0 } { $j < 3 } { incr j } {
				    # compute the magnitude of the resultant force in one node 
				    set total_bending [expr $total_bending + $mesh_bending($i,$j)*$mesh_bending($i,$j)]
				}
				set total_bending [expr sqrt($total_bending)]
				set kb_fmetric [expr $kb_fmetric + $total_bending]
			    }
			    set result $kb_fmetric
			}
		}
		if { $elastic_forces == "kal" || $elastic_forces == "total" || $elastic_forces == "kal-fmetric" } {
			set result 0

			# recover data of this object
			set object [lindex $oif_objects $objectID]
			set template_id [lindex $object 0]
			set template [lindex $oif_templates $template_id]
			set nnodes [lindex $template 0]
			set ntriangles [lindex $template 2]
			set kal [lindex $template 7]

			set start_id_of_nodes 0
			set start_id_of_triangles 0
			for {set i 0} {$i < $template_id} {incr i} {
			    set start_id_of_nodes [expr $start_id_of_nodes + [lindex [lindex $oif_templates $i] 0]]
			    set start_id_of_triangles [expr $start_id_of_triangles + [lindex [lindex $oif_templates $i] 2]]
			}

			# recover particles from this object's template  
			for {set i 0} {$i < $nnodes} {incr i} {
			    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_nodes+$i]]
			    set mesh_nodes($i,0) [lindex $node_triplet 0] 
			    set mesh_nodes($i,1) [lindex $node_triplet 1]
			    set mesh_nodes($i,2) [lindex $node_triplet 2]
			}

			# recover triangles from this object's template
			for {set i 0} {$i < $ntriangles} {incr i} {
			    set triangle [lindex $oif_template_triangles [expr $start_id_of_triangles+$i]]
			    set mesh_triangles($i,0) [lindex $triangle 0]
			    set mesh_triangles($i,1) [lindex $triangle 1]
			    set mesh_triangles($i,2) [lindex $triangle 2]
			}

			# initialize the mesh_local_area list
			for { set i 0 } { $i < $nnodes } { incr i } {
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_local_area($i,$j) 0.0
			    }
			}

			set start_id_of_particles [lindex $oif_object_starting_particles $objectID] 
			    
			for { set i 0 } { $i < $ntriangles } { incr i } {
			    # Take a triangle and copy the nodes of the triangle to pA, pB, pC
			    set pA $mesh_triangles($i,0)
			    set pB $mesh_triangles($i,1)
			    set pC $mesh_triangles($i,2)
			    
			    set origA [list $mesh_nodes($pA,0) $mesh_nodes($pA,1) $mesh_nodes($pA,2)]
			    set origB [list $mesh_nodes($pB,0) $mesh_nodes($pB,1) $mesh_nodes($pB,2)]
			    set origC [list $mesh_nodes($pC,0) $mesh_nodes($pC,1) $mesh_nodes($pC,2)]
			    
			    # Get the current position of the same three points
			    set currA [part [expr $start_id_of_particles + $pA] print pos]
			    set currAx [lindex $currA 0]
			    set currAy [lindex $currA 1]
			    set currAz [lindex $currA 2]
			    set currB [part [expr $start_id_of_particles + $pB] print pos]
			    set currBx [lindex $currB 0]
			    set currBy [lindex $currB 1]
			    set currBz [lindex $currB 2]
			    set currC [part [expr $start_id_of_particles + $pC] print pos]
			    set currCx [lindex $currC 0]
			    set currCy [lindex $currC 1]
			    set currCz [lindex $currC 2]
			    
			    # Calculate local area force for these three points
			    set orig_area [area_triangle origA origB origC]
			    set curr_area [area_triangle currA currB currC]
			    
			    set fTemp [calc_local_area_force $kal $currAx $currAy $currAz $currBx $currBy $currBz \
				       $currCx $currCy $currCz $orig_area $curr_area]
			    # first three numbers are force for point A, second three for point B, last three for point C
			    set f(0) [lindex $fTemp 0]
			    set f(1) [lindex $fTemp 1]
			    set f(2) [lindex $fTemp 2]
			    set f(3) [lindex $fTemp 3]
			    set f(4) [lindex $fTemp 4]
			    set f(5) [lindex $fTemp 5]
			    set f(6) [lindex $fTemp 6]
			    set f(7) [lindex $fTemp 7]
			    set f(8) [lindex $fTemp 8]
			      
			    # Save this local area force into the array of local area forces 
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_local_area($pA,$j) [expr $mesh_local_area($pA,$j) + $f([expr 0+$j])]
				set mesh_local_area($pB,$j) [expr $mesh_local_area($pB,$j) + $f([expr 3+$j])]
				set mesh_local_area($pC,$j) [expr $mesh_local_area($pC,$j) + $f([expr 6+$j])]
			    }
			    
			}
			# output here only if elastic forces are set to kal
			if { $elastic_forces == "kal" } {
			    # determine whether .vtk or other file type will be written
			    if { [ string first ".vtk" $output_file ] > -1 } {
			        # write local area forces into a .vtk file for visualization
				set out_f [open $output_file "w"]
			        puts $out_f "# vtk DataFile Version 3.0"
			        puts $out_f "Data"
			        puts $out_f "ASCII"
			        puts $out_f "DATASET POLYDATA"
			        puts $out_f "POINTS $nnodes float"
			        # writing points
			        for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
					puts $out_f "[part $i print pos]"
			        }
			        set ntriangles [lindex $template 2]
			        puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			        # extracts the starting position of the triangles for the template_id
			        set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			        # writing triangles
			        for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
					puts $out_f "3 [lindex $oif_template_triangles $i]"
			        }
			        # writing local area forces
			        puts $out_f "POINT_DATA $nnodes"
			        puts $out_f "SCALARS sample_scalars float 1"
			        puts $out_f "LOOKUP_TABLE default"
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_local_area 0.0
			        	for { set j 0 } { $j < 3 } { incr j } {
			        		# compute the magnitude of the resultant force to be written
						set total_local_area [expr $total_local_area + $mesh_local_area($i,$j)*$mesh_local_area($i,$j)]
					}
					set total_local_area [expr sqrt($total_local_area)]
					puts $out_f "$total_local_area"
			        }
			        close $out_f
			    } else {
			        # write local area forces into a text file
			        set out_f [open $output_file "w"]
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_local_area 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_local_area [expr $total_local_area + $mesh_local_area($i,$j)*$mesh_local_area($i,$j)]
					}
					set total_local_area [expr sqrt($total_local_area)]
					puts $out_f "$total_local_area"
			        }
			        close $out_f
			    }
			}
			if { $elastic_forces == "kal-fmetric" } {
			    set kal_fmetric 0.0 
			    for { set i 0 } { $i < $nnodes } { incr i } {
				set total_local_area 0.0
				for { set j 0 } { $j < 3 } { incr j } {
				    # compute the magnitude of the resultant force in one node 
				    set total_local_area [expr $total_local_area + $mesh_local_area($i,$j)*$mesh_local_area($i,$j)]
				}
				set total_local_area [expr sqrt($total_local_area)]
				set kal_fmetric [expr $kal_fmetric + $total_local_area]
			    }
			    set result $kal_fmetric
			}
		}
		if { $elastic_forces == "kag" || $elastic_forces == "total" || $elastic_forces == "kag-fmetric" } {
			set result 0

			# recover data of this object
			set object [lindex $oif_objects $objectID]
			set template_id [lindex $object 0]
			set template [lindex $oif_templates $template_id]
			set nnodes [lindex $template 0]
			set ntriangles [lindex $template 2]
			set kag [lindex $template 9]

			set start_id_of_nodes 0
			set start_id_of_triangles 0
			for {set i 0} {$i < $template_id} {incr i} {
			    set start_id_of_nodes [expr $start_id_of_nodes + [lindex [lindex $oif_templates $i] 0]]
			    set start_id_of_triangles [expr $start_id_of_triangles + [lindex [lindex $oif_templates $i] 2]]
			}

			# recover particles from this object's template  
			for {set i 0} {$i < $nnodes} {incr i} {
			    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_nodes+$i]]
			    set mesh_nodes($i,0) [lindex $node_triplet 0] 
			    set mesh_nodes($i,1) [lindex $node_triplet 1]
			    set mesh_nodes($i,2) [lindex $node_triplet 2]
			}

			# recover triangles from this object's template
			for {set i 0} {$i < $ntriangles} {incr i} {
			    set triangle [lindex $oif_template_triangles [expr $start_id_of_triangles+$i]]
			    set mesh_triangles($i,0) [lindex $triangle 0]
			    set mesh_triangles($i,1) [lindex $triangle 1]
			    set mesh_triangles($i,2) [lindex $triangle 2]
			}

			# initialize the mesh_global_area list
			for { set i 0 } { $i < $nnodes } { incr i } {
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_global_area($i,$j) 0.0
			    }
			}
			    
			# get the current and original global area
			set curr_gl_area [oif_object_analyze object-id $objectID surface-area]
			set orig_gl_area [lindex $template 13]

			set start_id_of_particles [lindex $oif_object_starting_particles $objectID] 

			for { set i 0 } { $i < $ntriangles } { incr i } {
			    # Take a triangle and copy the nodes of the triangle to pA, pB, pC
			    set pA $mesh_triangles($i,0)
			    set pB $mesh_triangles($i,1)
			    set pC $mesh_triangles($i,2)
			    
			    set origA [list $mesh_nodes($pA,0) $mesh_nodes($pA,1) $mesh_nodes($pA,2)]
			    set origB [list $mesh_nodes($pB,0) $mesh_nodes($pB,1) $mesh_nodes($pB,2)]
			    set origC [list $mesh_nodes($pC,0) $mesh_nodes($pC,1) $mesh_nodes($pC,2)]
			    
			    # Get the current position of the same three points
			    set currA [part [expr $start_id_of_particles + $pA] print pos]
			    set currAx [lindex $currA 0]
			    set currAy [lindex $currA 1]
			    set currAz [lindex $currA 2]
			    set currB [part [expr $start_id_of_particles + $pB] print pos]
			    set currBx [lindex $currB 0]
			    set currBy [lindex $currB 1]
			    set currBz [lindex $currB 2]
			    set currC [part [expr $start_id_of_particles + $pC] print pos]
			    set currCx [lindex $currC 0]
			    set currCy [lindex $currC 1]
			    set currCz [lindex $currC 2]
			    
			    # Calculate global area force for these three points
			    set fTemp [calc_global_area_force $kag $currAx $currAy $currAz $currBx $currBy $currBz \
				       $currCx $currCy $currCz $orig_gl_area $curr_gl_area]
			    # first three numbers are force for point A, second three for point B, last three for point C
			    set f(0) [lindex $fTemp 0]
			    set f(1) [lindex $fTemp 1]
			    set f(2) [lindex $fTemp 2]
			    set f(3) [lindex $fTemp 3]
			    set f(4) [lindex $fTemp 4]
			    set f(5) [lindex $fTemp 5]
			    set f(6) [lindex $fTemp 6]
			    set f(7) [lindex $fTemp 7]
			    set f(8) [lindex $fTemp 8]
			      
			    # Save this global area force into the array of global area forces 
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_global_area($pA,$j) [expr $mesh_global_area($pA,$j) + $f([expr 0+$j])]
				set mesh_global_area($pB,$j) [expr $mesh_global_area($pB,$j) + $f([expr 3+$j])]
				set mesh_global_area($pC,$j) [expr $mesh_global_area($pC,$j) + $f([expr 6+$j])]
			    }
			    
			}
			# output here only if elastic forces are set to kag
			if { $elastic_forces == "kag" } {
			    # determine whether .vtk or other file type will be written
			    if { [ string first ".vtk" $output_file ] > -1 } {
			        # write global area forces into a .vtk file for visualization
			        set out_f [open $output_file "w"]
			        puts $out_f "# vtk DataFile Version 3.0"
			        puts $out_f "Data"
			        puts $out_f "ASCII"
			        puts $out_f "DATASET POLYDATA"
			        puts $out_f "POINTS $nnodes float"
			        # writing points
			        for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
					puts $out_f "[part $i print pos]"
			        }
			        set ntriangles [lindex $template 2]
			        puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			        # extracts the starting position of the triangles for the template_id
			        set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			        # writing triangles
			        for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
					puts $out_f "3 [lindex $oif_template_triangles $i]"
			        }
			        # writing global area forces
			        puts $out_f "POINT_DATA $nnodes"
			        puts $out_f "SCALARS sample_scalars float 1"
			        puts $out_f "LOOKUP_TABLE default"
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_global_area 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_global_area [expr $total_global_area + $mesh_global_area($i,$j)*$mesh_global_area($i,$j)]
					}
					set total_global_area [expr sqrt($total_global_area)]
					puts $out_f "$total_global_area"
			        }
			        close $out_f
			    } else {
			        # write global area forces into a text file
			        set out_f [open $output_file "w"]
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_global_area 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_global_area [expr $total_global_area + $mesh_global_area($i,$j)*$mesh_global_area($i,$j)]
					}
					set total_global_area [expr sqrt($total_global_area)]
					puts $out_f "$total_global_area"
			        }
			        close $out_f
			    }
			}
			if { $elastic_forces == "kag-fmetric" } {
			    set kag_fmetric 0.0 
			    for { set i 0 } { $i < $nnodes } { incr i } {
				set total_global_area 0.0
				for { set j 0 } { $j < 3 } { incr j } {
				    # compute the magnitude of the resultant force in one node 
				    set total_global_area [expr $total_global_area + $mesh_global_area($i,$j)*$mesh_global_area($i,$j)]
				}
				set total_global_area [expr sqrt($total_global_area)]
				set kag_fmetric [expr $kag_fmetric + $total_global_area]
			    }
			    set result $kag_fmetric
			}
		}
		if { $elastic_forces == "kv" || $elastic_forces == "total" || $elastic_forces == "kv-fmetric" } {
			set result 0

			# recover data of this object
			set object [lindex $oif_objects $objectID]
			set template_id [lindex $object 0]
			set template [lindex $oif_templates $template_id]
			set nnodes [lindex $template 0]
			set ntriangles [lindex $template 2]
			set kv [lindex $template 11]

			set start_id_of_nodes 0
			set start_id_of_triangles 0
			for {set i 0} {$i < $template_id} {incr i} {
			    set start_id_of_nodes [expr $start_id_of_nodes + [lindex [lindex $oif_templates $i] 0]]
			    set start_id_of_triangles [expr $start_id_of_triangles + [lindex [lindex $oif_templates $i] 2]]
			}

			# recover particles from this object's template  
			for {set i 0} {$i < $nnodes} {incr i} {
			    set node_triplet [lindex $oif_template_nodes [expr $start_id_of_nodes+$i]]
			    set mesh_nodes($i,0) [lindex $node_triplet 0] 
			    set mesh_nodes($i,1) [lindex $node_triplet 1]
			    set mesh_nodes($i,2) [lindex $node_triplet 2]
			}

			# recover triangles from this object's template
			for {set i 0} {$i < $ntriangles} {incr i} {
			    set triangle [lindex $oif_template_triangles [expr $start_id_of_triangles+$i]]
			    set mesh_triangles($i,0) [lindex $triangle 0]
			    set mesh_triangles($i,1) [lindex $triangle 1]
			    set mesh_triangles($i,2) [lindex $triangle 2]
			}

			# initialize the mesh_volume list
			for { set i 0 } { $i < $nnodes } { incr i } {
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_volume($i,$j) 0.0
			    }
			}
			    
			# compute the current and original volume
			set curr_volume [oif_object_analyze object-id $objectID volume]
			set orig_volume [lindex $template 14]

			set start_id_of_particles [lindex $oif_object_starting_particles $objectID] 

			for { set i 0 } { $i < $ntriangles } { incr i } {
			    # Take a triangle and copy the nodes of the triangle to pA, pB, pC
			    set pA $mesh_triangles($i,0)
			    set pB $mesh_triangles($i,1)
			    set pC $mesh_triangles($i,2)
			    
			    set origA [list $mesh_nodes($pA,0) $mesh_nodes($pA,1) $mesh_nodes($pA,2)]
			    set origB [list $mesh_nodes($pB,0) $mesh_nodes($pB,1) $mesh_nodes($pB,2)]
			    set origC [list $mesh_nodes($pC,0) $mesh_nodes($pC,1) $mesh_nodes($pC,2)]
			    
			    # Get the current position of the same three points
			    set currA [part [expr $start_id_of_particles + $pA] print pos]
			    set currAx [lindex $currA 0]
			    set currAy [lindex $currA 1]
			    set currAz [lindex $currA 2]
			    set currB [part [expr $start_id_of_particles + $pB] print pos]
			    set currBx [lindex $currB 0]
			    set currBy [lindex $currB 1]
			    set currBz [lindex $currB 2]
			    set currC [part [expr $start_id_of_particles + $pC] print pos]
			    set currCx [lindex $currC 0]
			    set currCy [lindex $currC 1]
			    set currCz [lindex $currC 2]
			    
			    # Calculate volume force for these three points
			    set fTemp [calc_volume_force $kv $currAx $currAy $currAz $currBx $currBy $currBz \
				       $currCx $currCy $currCz $orig_volume $curr_volume]
			    set f(0) [lindex $fTemp 0]
			    set f(1) [lindex $fTemp 1]
			    set f(2) [lindex $fTemp 2]
			    
			    # Save this global area force into the array of volume forces 
			    for { set j 0 } { $j < 3 } { incr j } {
				set mesh_volume($pA,$j) [expr $mesh_volume($pA,$j) + $f($j)]
				set mesh_volume($pB,$j) [expr $mesh_volume($pB,$j) + $f($j)]
				set mesh_volume($pC,$j) [expr $mesh_volume($pC,$j) + $f($j)]
			    }
			    
			}
			# output here only if elastic forces are set to kv
			if { $elastic_forces == "kv" } {
			    # determine whether .vtk or other file type will be written
			    if { [ string first ".vtk" $output_file ] > -1 } {
			        # write volume forces into a .vtk file for visualization
			        set out_f [open $output_file "w"]
			        puts $out_f "# vtk DataFile Version 3.0"
			        puts $out_f "Data"
			        puts $out_f "ASCII"
			        puts $out_f "DATASET POLYDATA"
			        puts $out_f "POINTS $nnodes float"
			        # writing points
			        for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
					puts $out_f "[part $i print pos]"
			        }
			        set ntriangles [lindex $template 2]
			        puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			        # extracts the starting position of the triangles for the template_id
			        set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			        # writing triangles
			        for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
					puts $out_f "3 [lindex $oif_template_triangles $i]"
			        }
			        # writing volume forces
			        puts $out_f "POINT_DATA $nnodes"
			        puts $out_f "SCALARS sample_scalars float 1"
			        puts $out_f "LOOKUP_TABLE default"
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_volume 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_volume [expr $total_volume + $mesh_volume($i,$j)*$mesh_volume($i,$j)]
					}
					set total_volume [expr sqrt($total_volume)]
					puts $out_f "$total_volume"
			        }
			        close $out_f
			    } else {
			        # write volume forces into a text file
			        set out_f [open $output_file "w"]
			        for { set i 0 } { $i < $nnodes } { incr i } {
					set total_volume 0.0
					for { set j 0 } { $j < 3 } { incr j } {
						# compute the magnitude of the resultant force to be written
						set total_volume [expr $total_volume + $mesh_volume($i,$j)*$mesh_volume($i,$j)]
					}
					set total_volume [expr sqrt($total_volume)]
					puts $out_f "$total_volume"
			        }
			        close $out_f
			    }
			}
			if { $elastic_forces == "kv-fmetric" } {
			    set kv_fmetric 0.0 
			    for { set i 0 } { $i < $nnodes } { incr i } {
				set total_volume 0.0
				for { set j 0 } { $j < 3 } { incr j } {
				    # compute the magnitude of the resultant force in one node 
				    set total_volume [expr $total_volume + $mesh_volume($i,$j)*$mesh_volume($i,$j)]
				}
				set total_volume [expr sqrt($total_volume)]
				set kv_fmetric [expr $kv_fmetric + $total_volume]
			    }
			    set result $kv_fmetric
			}
		}
		if { $elastic_forces == "total" } {
		    # everything has already been computed, it will be just put together here
		    # determine whether .vtk or other file type will be written
		    if { [ string first ".vtk" $output_file ] > -1 } {
			# write total elastic forces into a .vtk file for visualization
			set out_f [open $output_file "w"]
			puts $out_f "# vtk DataFile Version 3.0"
			puts $out_f "Data"
			puts $out_f "ASCII"
			puts $out_f "DATASET POLYDATA"
			puts $out_f "POINTS $nnodes float"
			# writing points
			for { set i $start_id_of_particles } { $i < [expr $start_id_of_particles + $nnodes] } { incr i } {
			    puts $out_f "[part $i print pos]"
			}
			set ntriangles [lindex $template 2]
			puts $out_f "TRIANGLE_STRIPS $ntriangles [expr 4*$ntriangles]"
			# extracts the starting position of the triangles for the template_id
			set start_id_of_triangles [lindex $oif_template_starting_triangles $template_id]
			# writing triangles
			for { set i $start_id_of_triangles } { $i < [expr $start_id_of_triangles + $ntriangles] } { incr i } {
			    puts $out_f "3 [lindex $oif_template_triangles $i]"
			}
			# writing total elastic forces
			puts $out_f "POINT_DATA $nnodes"
			puts $out_f "SCALARS sample_scalars float 1"
			puts $out_f "LOOKUP_TABLE default"
			for { set i 0 } { $i < $nnodes } { incr i } {
			    # adding up all alastic forces in each node
			    for { set j 0 } { $j < 3 } { incr j } {
				set el($j) [expr $mesh_stretching($i,$j) + $mesh_bending($i,$j) + $mesh_local_area($i,$j) \
					    + $mesh_global_area($i,$j) + $mesh_volume($i,$j)]
			    }
			    # compute the magnitude of the resultant force to be written
			    set total_elasticity [expr sqrt($el(0)*$el(0)+$el(1)*$el(1)+$el(2)*$el(2))]
			    puts $out_f "$total_elasticity"
			}
			close $out_f
		    } else {
			# write elastic forces into a text file
			set out_f [open $output_file "w"]
			for { set i 0 } { $i < $nnodes } { incr i } {
			    # adding up all alastic forces in each node
			    for { set j 0 } { $j < 3 } { incr j } {
				set el($j) [expr $mesh_stretching($i,$j) + $mesh_bending($i,$j) + $mesh_local_area($i,$j) \
					    + $mesh_global_area($i,$j) + $mesh_volume($i,$j)]
			    }
			    # compute the magnitude of the resultant force to be written
			    set total_elasticity [expr sqrt($el(0)*$el(0)+$el(1)*$el(1)+$el(2)*$el(2))]
			    puts $out_f "$total_elasticity"
			}
			close $out_f
		    }
		    set result 1
		}

		if { $result == -1 } {
			puts "Argument of elastic-forces computation must be one of these: ks, kb, kal, kag, kv, total, \
			ks-fmetric, kb-fmetric, kal-fmetric, kag-fmetric, kv-fmetric"
		}
		return $result
	} 
}

proc oif_mesh_analyze { args } {
	# acces global variables defined in oif_init:
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
	set repair 0
	set corr_file ""
	set method -1
	set flip -1
	set shift_node_ids -1

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments at the end of this section
    set pos 0
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"nodes-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_mesh_analyze: nodes-file missing"
					break
				}
				set mesh_nodes_file [lindex $args $pos]
				incr pos
			}
			"triangles-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_mesh_analyze: triangles-file missing"
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
					puts "error in oif_mesh_analyze: output file for repaired triangulation missing"
					break
				}
				set repair 1
				set corr_file [lindex $args $pos]
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_mesh_analyze: method id for repair missing"
					break
				}
				set method [lindex $args $pos]
				incr pos
			}
			"flip" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_mesh_analyze: output file for flipped triangulation missing"
					break
				}
				set corr_file [lindex $args $pos]
				set flip 1
				incr pos
			}
			"shift-node-ids" {
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_mesh_analyze: output file for shifted nodes missing"
					break
				}
				set corr_file [lindex $args $pos]
				set shift_node_ids 1
				incr pos
			}
			default { 
				puts "error in oif_mesh_analyze: incorrect keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $mesh_nodes_file == "" } { set mandatory 0 }
	if { $mesh_triangles_file == "" } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_mesh_analyze: mandatory argument(s) missing" 
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
				puts "wrong orientation of triangle at line $i"
			}
			set P1 [lreplace $P1 0 2]
			set P2 [lreplace $P2 0 2]
			set P3 [lreplace $P3 0 2]

		}

		# Method 2: Second check controls whether all couples of triangles with the same edge have the same orientation. Not implementedn yet

	}

	if { $repair == 1} {
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
			# set tmp_res [expr $origin_X*$norm0 + $origin_Y*$norm1 + $origin_Z*$norm2 + $tmp_d]
			# puts "tmp_d: $tmp_d"
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
		# Method 2: Second check controls whether all couples of triangles with the same edge have the same orientation. Not implemented yet
		# TO BE IMPLEMENTED
	    }
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
	
	if { $shift_node_ids == 1} {
		set fcorr [open $corr_file w]
		set fp [open $mesh_triangles_file r]
		set file_data [read $fp]
		close $fp
		set data [split $file_data "\n"]
		foreach line $data {
			if { [llength $line] == 3 } {
				set v0 [expr [lindex $line 0]-1]
				set v1 [expr [lindex $line 1]-1]
				set v2 [expr [lindex $line 2]-1]
				puts $fcorr "$v0 $v1 $v2"
			}
		}
		close $fcorr
	}

	
}

proc oif_obj_obj_interactions { args } {
	# access global variables
	global oif_n_objects
	global oif_nparticles
	global oif_ntriangles
	global oif_nedges
	global oif_firstBondId
	global oif_firstPartId
	global oif_objects
	global oif_object_starting_particles
	global oif_template_starting_triangles
	global oif_object_object_interactions
	set n_args 0
		# counts the number of arguments
	foreach arg $args {
		incr n_args
    }
	if { $n_args == 0 } {
		puts "oif_obj_obj_interactions: Mandatory arguments are: cut-off value" 
		return 0
	}

	##### reading the arguments. some of them are mandatory. we check for the mandatory arguments ad the end of this section
    set pos 0
    set cut_off -1
    set output_file ""
    while { $pos < $n_args } {
		switch [lindex $args $pos] {
			"cut-off" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_obj_obj_interactions: missing cut-off value"
					break
				}
				set cut_off [lindex $args $pos]
				incr pos
			}
			"output-file" {  
				incr pos
				if { $pos >= $n_args } { 
					puts "error in oif_obj_obj_interactions: missing output filename"
					break
				}
				set output_file [lindex $args $pos]
				incr pos
			}
			default { 
				puts "error in oif_obj_obj_interactions: incorrect keyword" 
				set pos $n_args
			}
		}  
	}

# checking whether all mandatory arguments have been given
	set mandatory 1
	if { $cut_off == -1 } { set mandatory 0 }

	if { $mandatory == 0 } { 
		puts "error in oif_obj_obj_interactions: mandatory argument(s) missing" 
		return
	}
# gets coordinates of all objects
	for { set obj 0 } { $obj < $oif_n_objects } {incr obj } {
		set orig [ oif_object_analyze object-id $obj approx-pos]
		set origX($obj) [lindex $orig 0]
		set origY($obj) [lindex $orig 1]
		set origZ($obj) [lindex $orig 2]
		
	}
#set pairwise distances between all objects
	puts "setting pairwise distances between objects"
	set min_dist 1000000
	for { set objA 0 } { $objA < [expr $oif_n_objects ] } {incr objA } {
		for { set objB [expr $objA + 1] } { $objB < [expr $oif_n_objects] } {incr objB } {
			set dx [expr $origX($objA) - $origX($objB)]
			set dy [expr $origY($objA) - $origY($objB)]
			set dz [expr $origZ($objA) - $origZ($objB)]
			set mut_dist [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]
			set obj_mutual_distances($objA,$objB) $mut_dist
			if { $mut_dist < $min_dist } { set min_dist $mut_dist } 
		}
	}
	
	if {$output_file != "" } {
		puts "saving pairwise distances into file"
		set fout [open $output_file w]
		for { set objA 0 } { $objA < [expr $oif_n_objects ] } {incr objA } {
			for { set objB [expr $objA + 1] } { $objB < [expr $oif_n_objects ] } {incr objB } {
				puts $fout "$objA - $objB $obj_mutual_distances($objA,$objB)"
			}
		}
		close $fout
	}
	puts "Before"
	puts "$oif_object_object_interactions"
	
	for { set objA 0 } { $objA < [expr $oif_n_objects ] } {incr objA } {
		if {$output_file != "" } {
			set fout [open "state-$output_file" a]
			puts $fout "checking $objA out of $oif_n_objects"
			close $fout
		}
		for { set objB [expr $objA + 1] } { $objB < [expr $oif_n_objects] } {incr objB } {
#			puts "[lsearch $oif_object_object_interactions { $objA $objB }]"
#			puts "$obj_mutual_distances($objA,$objB), $cut_off"
			if { $obj_mutual_distances($objA,$objB) < $cut_off } {
				if { [lsearch $oif_object_object_interactions [list $objA $objB]] == -1 } {
#					puts "[lsearch $oif_object_object_interactions [list $objA $objB]]"
					lappend oif_object_object_interactions [list $objA $objB]
#					puts "[lsearch $oif_object_object_interactions [list $objA $objB]]"
					#add interaction
					inter $objA $objB soft-sphere 0.01 1.2 0.2 0.0
					#puts "Adding soft sphere $objA $objB"
				}
			}
		}
	}
	
	puts "After:"
	puts "$oif_object_object_interactions"

	foreach line $oif_object_object_interactions {
		incr n_obj_inter([lindex $line 0])
		incr n_obj_inter([lindex $line 1])
    }
	set n_obj_obj_inter 0
	if {$output_file != "" } {
		set fout [open "list-$output_file" a]

		puts $fout "Number of non-bonded interactions for each object"
		if { [llength $oif_object_object_interactions] > 0 } {
			for { set obj 0 } { $obj < $oif_n_objects } { incr obj } {
				puts $fout "$n_obj_inter($obj)"
				set n_obj_obj_inter [expr $n_obj_obj_inter + $n_obj_inter($obj)]
			}	
			set n_obj_obj_inter [expr $n_obj_obj_inter/2]
		}
		puts $fout "Number of all pairs of non-bonded interactions: $n_obj_obj_inter"
		puts $fout "All possible pairwise interactions: [expr $oif_n_objects * ($oif_n_objects - 1)/2]"
		close $fout
	}
}
