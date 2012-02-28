# Copyright (C) 2012 The ESPResSo project
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

proc dielectric_wall { args } {
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  if { ![ info exists n_induced_charges ] } {
    set n_induced_charges 0 
    set icc_areas [ list ]
    set icc_normals [ list ] 
    set icc_epsilons [ list ] 
    set icc_sigmas [ list ]
  }
  set dist 0
  set nx 0
  set ny 0
  set nz 0
  set res 0
  set sigma 0
  set type 0 
  set eps 0.
  for { set argno 0 } { $argno < [ llength $args ] } { incr argo } {
    if { [ lindex $args $argno ] == "dist" } {
      incr argno
      set dist [ expr 1.0 * [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "normal" } {
      incr argno
      set nx [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      set ny [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      set nz [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "res" } {
      incr argno
      set res [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "sigma" } {
      incr argno
      set sigma [ expr 1.0* [ lindex $args $argno ]] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "eps" } {
      incr argno
      set eps [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "type" } {
      incr argno
      set type [ expr 1.0* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    error "did not understand arg [ lindex args $argno ]"
  }
  set box_l_x [ lindex [ setmd box_l ] 0 ]
  set box_l_y [ lindex [ setmd box_l ] 1 ]
  set box_l_z [ lindex [ setmd box_l ] 2 ]
  set counter 0
# now construct two normal vector that are perpendicular to each other and to the orientation vector:
  if [ expr !(abs(abs($nx) - 1. )< 0.01) ] {
    set e_1_x 0
    set e_1_y [ expr - $nz ]
    set e_1_z $ny
  } elseif { abs((abs($ny) - 1) < 0.01) } {
    set e_1_x $nz
    set e_1_y 0
    set e_1_z [ expr - $nx ]
  } elseif { !(abs($nz) - 1 < 0.01) } {
    set e_1_x $ny
    set e_1_y [ expr -$nx ]
    set e_1_z 0
  } else {
    set e_1_x 0
    set e_1_y [ expr - $nz ]
    set e_1_z $ny
  }
  puts "e_1 is $e_1_x $e_1_y $e_1_z"
# and normalize
  set sq_sum [ expr sqrt($e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z) ] 
  set e_1_x [ expr $e_1_x/$sq_sum ]
  set e_1_y [ expr $e_1_y/$sq_sum ]
  set e_1_z [ expr $e_1_z/$sq_sum ]
  set sq_sum [ expr $e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z ] 
# cross product
  set e_2_x [ expr $ny*$e_1_z - $nz*$e_1_y ]
  set e_2_y [ expr $nz*$e_1_x - $nx*$e_1_z ]
  set e_2_z [ expr $nx*$e_1_y - $ny*$e_1_x ]
  set sq_sum [ expr $e_2_x*$e_2_x +  $e_2_y*$e_2_y + $e_2_z*$e_2_z ] 

# We now go from -box_l to box_l in e1 and e2 direction and put particles 
# in all positions that are in the within the box.
  set n_charges_l [ expr int(($box_l_x + $box_l_y + $box_l_z)/$res) ]
  set l_incr $res
  for { set i -$n_charges_l } { $i <= $n_charges_l } { incr i } {
    for { set j -$n_charges_l } { $j <= $n_charges_l } { incr j } {
      set l1 $i*$l_incr
      set l2 $j*$l_incr
      set posx [ expr $nx * $dist + $l1*$e_1_x + $l2*$e_2_x ]
      set posy [ expr $ny * $dist + $l1*$e_1_y + $l2*$e_2_y ]
      set posz [ expr $nz * $dist + $l1*$e_1_z + $l2*$e_2_z ]
      if { $posx >= 0 && $posx < $box_l_x &&  $posy >= 0 && $posy < $box_l_y &&  $posz >= 0 && $posz < $box_l_z } {
        part $n_induced_charges pos $posx $posy $posz q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] fix 1 1 1 type $type
        incr n_induced_charges
        lappend icc_normals [ list $nx $ny $nz ] 
        incr counter
      }
    }
  }
  for { set i 0 } { $i < $counter } { incr i } {
    lappend icc_areas [ expr $res*$res ]
    lappend icc_sigmas $sigma 
    lappend icc_epsilons $eps 
  }


}


proc dielectric_sphere { args } {
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  if { ![ info exists n_induced_charges ] } {
    set n_induced_charges 0 
    set icc_areas [ list ]
    set icc_normals [ list ] 
    set icc_epsilons [ list ] 
    set icc_sigmas [ list ]
  }
  set radius 0
  set cx 0
  set cy 0
  set cz 0
  set direction 1
  set res 0
  set sigma 0
  set type 0 
  for { set argno 0 } { $argno < [ llength $args ] } { incr argo } {
    if { [ lindex $args $argno ] == "radius" } {
      incr argno
      set radius [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "center" } {
      incr argno
      set cx [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      set cy [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      set cz [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "res" } {
      incr argno
      set res [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "direction" } {
      incr argno
      set direction [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "sigma" } {
      incr argno
      set sigma [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "eps" } {
      incr argno
      set eps [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "type" } {
      incr argno
      set type [ expr 1.* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    error "did not understand arg [ lindex args $argno ]"
  }
  set pi 3.1415
  set n_half_equator [ expr round($pi * $radius/$res) ]

  set ic_counter 0

  set theta 0.
  set incr_theta [ expr $pi/$n_half_equator ] 
  for { set i 0 } { $i <= $n_half_equator } { incr i } {
    set phi 0.
    set n_circle [ expr 2*round($n_half_equator*sin($theta)) ]
    set incr_phi [expr 2*$pi/$n_circle]
    for { set j 0 } { $j < $n_circle } { incr j} {
       set x [ expr  $radius*cos($phi)*sin($theta) ]
       set y [ expr  $radius*sin($phi)*sin($theta) ]
       set z [ expr  $radius*cos($theta) ]
       set phi [ expr $phi+$incr_phi ]
       part [ expr $n_induced_charges ] pos [expr $cx+$x] [expr $cy+$y] [expr $cy+$z] q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] type $type fix 1 1 1
       incr n_induced_charges
       incr ic_counter
       lappend icc_normals [ list  [ expr $direction*$x/$radius ] [ expr $direction* $y/$radius ] [ expr $direction* $z/$radius ] ]
    }
    set theta [ expr $theta + $incr_theta ]
  }
  for { set i 0 } { $i < $ic_counter } { incr i } {
    lappend icc_areas [ expr 4*$pi*$radius*$radius/double($ic_counter) ]
    lappend icc_epsilons $eps
    lappend icc_sigmas $sigma 
  }
}





proc dielectric_cylinder { args } {
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  if { ![ info exists n_induced_charges ] } {
    set n_induced_charges 0 
    set icc_areas [ list ]
    set icc_normals [ list ] 
    set icc_epsilons [ list ] 
    set icc_sigmas [ list ]
  }
  set center_x 0
  set center_y 0
  set center_z 0
  set axis_x 0
  set axis_y 0
  set axis_z 0
  set radius 0
  set length 0
  set n_length 0
  set n_circle 0
  set eps 0
  set covers 0
  set direction 1
  set sigma 0
  set type 0
  for { set argno 0 } { $argno < [ llength $args ] } { incr argo } {
    if { [ lindex $args $argno ] == "center" } {
      incr argno
      set center_x [ lindex $args $argno ] 
      incr argno
      set center_y [ lindex $args $argno ] 
      incr argno
      set center_z [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "axis" } {
      incr argno
      set axis_x [ lindex $args $argno ] 
      incr argno
      set axis_y [ lindex $args $argno ] 
      incr argno
      set axis_z [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "radius" } {
      incr argno
      set radius [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "length" } {
      incr argno
      set length [ expr 2* [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "res" } {
      incr argno
      set res [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "eps" } {
      incr argno
      set eps [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "type" } {
      incr argno
      set type [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "covers" } {
      incr argno
      set covers [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "direction" } {
      incr argno
      set direction [ lindex $args $argno ] 
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "sigma" } {
      incr argno
      set sigma [ lindex $args $argno ] 
      incr argno
      continue
    }
    error "did not understand argument [ lindex $args $argno ]"
  }

  set pi 3.1415
  set n_circle [ expr round( 2*$pi*$radius/$res) ]
  set n_length [ expr round( $length/$res) ]
  set particle_counter 0
# normalize orientation vector
  set sq_sum [ expr sqrt(1.0*$axis_x*$axis_x +  $axis_y*$axis_y + $axis_z*$axis_z) ] 
  set axis_x [ expr $axis_x/$sq_sum ]
  set axis_y [ expr $axis_y/$sq_sum ]
  set axis_z [ expr $axis_z/$sq_sum ]
# The polar angle is changed by incr phi 
  set incr_phi [ expr 2*$pi / $n_circle ]
# the cylinder is constructed by iteration along the axis, starting at start_* and then in circles
  set incr_length_x [ expr $axis_x * $length / $n_length]
  set incr_length_y [ expr $axis_y * $length / $n_length]
  set incr_length_z [ expr $axis_z * $length / $n_length]
  set start_x [ expr $center_x - $axis_x * $length/2. ]
  set start_y [ expr $center_y - $axis_y * $length/2. ]
  set start_z [ expr $center_z - $axis_z * $length/2. ]
  set circle_center_x $start_x
  set circle_center_y $start_y
  set circle_center_z $start_z
# now construct two normal vector that are perpendicular to each other and to the orientation vector:
  if [ expr !abs($axis_x - 1.) < 0.01 ] {
    set e_1_x 0
    set e_1_y [ expr - $axis_z ]
    set e_1_z $axis_y
  } elseif { !abs($axis_y - 1) < 0.01 } {
    set e_1_x $axis_z
    set e_1_y 0
    set e_1_z [ expr - $axis_x ]
  } elseif { !abs($axis_z - 1) < 0.01 } {
    set e_1_x $axis_y
    set e_1_y [ expr -$axis_x ]
    set e_1_z 0
  } else {
    set e_1_x 0
    set e_1_y [ expr - $axis_z ]
    set e_1_z $axis_y
  }
# and normalize
  set sq_sum [ expr sqrt($e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z) ] 
  set e_1_x [ expr $e_1_x/$sq_sum ]
  set e_1_y [ expr $e_1_y/$sq_sum ]
  set e_1_z [ expr $e_1_z/$sq_sum ]
  set sq_sum [ expr $e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z ] 
# cross product
  set e_2_x [ expr $axis_y*$e_1_z - $axis_z*$e_1_y ]
  set e_2_y [ expr $axis_z*$e_1_x - $axis_x*$e_1_z ]
  set e_2_z [ expr $axis_x*$e_1_y - $axis_y*$e_1_x ]
  set sq_sum [ expr $e_2_x*$e_2_x +  $e_2_y*$e_2_y + $e_2_z*$e_2_z ] 
# first the curved cylinder wall
  set phi 0.
  for { set i 0 } { $i <= $n_length} { incr i} {
    for { set j 0 } { $j < $n_circle } { incr j} {
      set x [ expr $circle_center_x + $radius*cos($phi)*$e_1_x + $radius*sin($phi)*$e_2_x ]
      set y [ expr $circle_center_y + $radius*cos($phi)*$e_1_y + $radius*sin($phi)*$e_2_y ]
      set z [ expr $circle_center_z + $radius*cos($phi)*$e_1_z + $radius*sin($phi)*$e_2_z ]
      part [ expr $n_induced_charges ] pos $x $y $z type $type q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] type $type
      set nx [ expr $direction*(cos($phi)*$e_1_x + sin($phi)*$e_2_x) ]
      set ny [ expr $direction*(cos($phi)*$e_1_y + sin($phi)*$e_2_y) ]
      set nz [ expr $direction*(cos($phi)*$e_1_z + sin($phi)*$e_2_z) ]
      lappend icc_normals [ list $nx $ny $nz ]
      incr particle_counter
      incr n_induced_charges
      set phi $phi+$incr_phi
    }
    set circle_center_x [ expr $circle_center_x + $incr_length_x ]
    set circle_center_y [ expr $circle_center_y + $incr_length_y ]
    set circle_center_z [ expr $circle_center_z + $incr_length_z ]
  }
# now the covers
  if {$covers} {
    set n_r [ expr floor($radius/$res) ]
    set r_incr [ expr $radius / $n_r ]
    set circle_center_x $start_x
    set circle_center_y $start_y
    set circle_center_z $start_z
    set this_radius 0
    for { set j 0 } { $j < $n_r } { incr j } { 
      set n_circle [ expr round( 2*$pi*$this_radius/$res) ]
      set incr_phi [ expr 2*$pi / $n_circle ]
      for { set i 0 } { $i < $n_circle } { incr i } {
        set x [ expr $circle_center_x + $this_radius*cos($phi)*$e_1_x + $this_radius*sin($phi)*$e_2_x ]
        set y [ expr $circle_center_y + $this_radius*cos($phi)*$e_1_y + $this_radius*sin($phi)*$e_2_y ]
        set z [ expr $circle_center_z + $this_radius*cos($phi)*$e_1_z + $this_radius*sin($phi)*$e_2_z ]
        part [ expr $n_induced_charges ] pos $x $y $z type $type q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ]
        set nx [ expr -$direction*$axis_x  ]
        set ny [ expr -$direction*$axis_y  ]
        set nz [ expr -$direction*$axis_z  ]
        lappend icc_normals [ list $nx $ny $nz ] 
        incr particle_counter
        incr n_induced_charges
        set phi $phi+$incr_phi
      }
      set this_radius [ expr $this_radius + $r_incr ] 
    }
    set circle_center_x [ expr $center_x + $axis_x * $length/2. ]
    set circle_center_y [ expr $center_y + $axis_y * $length/2. ]
    set circle_center_z [ expr $center_z + $axis_z * $length/2. ]
    set this_radius 0
    for { set j 0 } { $j < $n_r } { incr j } { 
      set n_circle [ expr round( 2*$pi*$this_radius/$res) ]
      set incr_phi [ expr 2*$pi / $n_circle ]
      for { set i 0 } { $i < $n_circle } { incr i } {
        set x [ expr $circle_center_x + $this_radius*cos($phi)*$e_1_x + $this_radius*sin($phi)*$e_2_x ]
        set y [ expr $circle_center_y + $this_radius*cos($phi)*$e_1_y + $this_radius*sin($phi)*$e_2_y ]
        set z [ expr $circle_center_z + $this_radius*cos($phi)*$e_1_z + $this_radius*sin($phi)*$e_2_z ]
        part [ expr $n_induced_charges ] pos $x $y $z type $type q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] 
        set nx [ expr +$direction*$axis_x  ]
        set ny [ expr +$direction*$axis_y  ]
        set nz [ expr +$direction*$axis_z  ]
        lappend icc_normals [ list $nx $ny $nz ] 
        incr particle_counter
        incr n_induced_charges
        set phi $phi+$incr_phi
      }
      set this_radius [ expr $this_radius + $r_incr ] 
    }
  } 
  for { set i 0 } { $i < $particle_counter } { incr i } { 
    if { $covers } {
      lappend icc_areas [ expr ( 2*$pi*$radius*$length + $pi*$radius*$radius )/$particle_counter ]
    } else { 
      lappend icc_areas [ expr ( 2*$pi*$radius*$length )/$particle_counter ]
    }
    lappend icc_epsilons $eps
    lappend icc_sigmas $sigma 
  }
}







proc sign x { 
  if { $x >= 0 } { 
    return 1
  } else {
    return -1
  }
}


proc dielectric_pore { args } { 
  global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
  if { ![ info exists n_induced_charges ] } {
    set n_induced_charges 0 
    set icc_areas [ list ]
    set icc_normals [ list ] 
    set icc_epsilons [ list ] 
    set icc_sigmas [ list ]
  }
  set smoothing_radius 1
  set r1 0.
  set r2 0.
  set length 0.
  set res 0.
  set type 0
  set eps 0.
  set sigma 0.
  set box_l_x [ lindex [ setmd box_l ] 0 ]
  set box_l_y [ lindex [ setmd box_l ] 1 ]
  set box_l_z [ lindex [ setmd box_l ] 2 ]
  set max_radius [ expr sqrt($box_l_x*$box_l_x+ $box_l_y*$box_l_y + $box_l_z*$box_l_z) ]

  for { set argno 0 } { $argno < [ llength $args ] } { incr argo } {
    if { [ lindex $args $argno ] == "center" } {
      incr argno
      set center_x [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      set center_y [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      set center_z [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "axis" } {
      incr argno
      set axis_x [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      set axis_y [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      set axis_z [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "radius" } {
      incr argno
      set r1 [ expr 1.0*[ lindex $args $argno ] ]
      set r2 $r1
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "radii" } {
      incr argno
      set r1 [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      set r2 [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "length" } {
      incr argno
      set length [ expr 2.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "res" } {
      incr argno
      set res [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "eps" } {
      incr argno
      set eps [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "type" } {
      incr argno
      set type [ expr [ lindex $args $argno ] ]
      incr argno
      continue
    }
    if { [ lindex $args $argno ] == "smoothing_radius" } {
      incr argno
      set smoothing_radius [ expr 1.0*[ lindex $args $argno ] ]
      incr argno
      continue
    }
    error "did not understand arg [ lindex args $argno ]"
  }
  set pi 3.1415
  set particle_counter 0
# normalize orientation vector
  set sq_sum [ expr sqrt(1.0*$axis_x*$axis_x +  $axis_y*$axis_y + $axis_z*$axis_z) ] 
  set axis_x [ expr $axis_x/$sq_sum ]
  set axis_y [ expr $axis_y/$sq_sum ]
  set axis_z [ expr $axis_z/$sq_sum ]
# The polar angle is changed by incr phi 
# the cylinder is constructed by iteration along the axis, starting at start_* and then in circles
  set start_x [ expr $center_x - $axis_x * $length/2. ]
  set start_y [ expr $center_y - $axis_y * $length/2. ]
  set start_z [ expr $center_z - $axis_z * $length/2. ]
  set circle_center_x $start_x
  set circle_center_y $start_y
  set circle_center_z $start_z
# now construct two normal vector that are perpendicular to each other and to the orientation vector:
  if [ expr !abs($axis_x - 1.) < 0.01 ] {
    set e_1_x 0
    set e_1_y [ expr - $axis_z ]
    set e_1_z $axis_y
  } elseif { !abs($axis_y - 1) < 0.01 } {
    set e_1_x $axis_z
    set e_1_y 0
    set e_1_z [ expr - $axis_x ]
  } elseif { !abs($axis_z - 1) < 0.01 } {
    set e_1_x $axis_y
    set e_1_y [ expr -$axis_x ]
    set e_1_z 0
  } else {
    set e_1_x 0
    set e_1_y [ expr - $axis_z ]
    set e_1_z $axis_y
  }
# and normalize
  set sq_sum [ expr sqrt($e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z) ] 
  set e_1_x [ expr $e_1_x/$sq_sum ]
  set e_1_y [ expr $e_1_y/$sq_sum ]
  set e_1_z [ expr $e_1_z/$sq_sum ]
  set sq_sum [ expr $e_1_x*$e_1_x +  $e_1_y*$e_1_y + $e_1_z*$e_1_z ] 
# cross product
  set e_2_x [ expr $axis_y*$e_1_z - $axis_z*$e_1_y ]
  set e_2_y [ expr $axis_z*$e_1_x - $axis_x*$e_1_z ]
  set e_2_z [ expr $axis_x*$e_1_y - $axis_y*$e_1_x ]
  set sq_sum [ expr $e_2_x*$e_2_x +  $e_2_y*$e_2_y + $e_2_z*$e_2_z ] 
  set slope [ expr ($r2 - $r1)/$length ]
  set rav [ expr ($r2 + $r1)/2 ]
  set p1 -1e5
  set p2 +1e5
  set length [ expr $length /2 ]

  set c1_z [ expr  - ($length - $smoothing_radius) ]
  set c2_z [ expr    ($length - $smoothing_radius)]

  set c1_r [ expr $r1 + $smoothing_radius ]
  set c2_r [ expr $r2 + $smoothing_radius ]
  set slope [ expr ($c2_r-$c1_r)/($c2_z-$c1_z) ]
  set sina [ expr ($r2-$r1)/sqrt( pow(2*($length - $smoothing_radius),2) + pow($r1-$r2,2)) ]
  set cosa [ expr sqrt(1-$sina*$sina) ]
  
#  set z_left  [ expr $c1_z + [ sign $slope ] * sqrt($slope*$slope/(1+$slope*$slope))*$smoothing_radius]
#  set z_right [ expr $c2_z + [ sign $slope ] * sqrt($slope*$slope/(1+$slope*$slope))*$smoothing_radius]
  set z_left  [ expr $c1_z + $sina * $smoothing_radius ]
  set z_right [ expr $c2_z + $sina * $smoothing_radius ]

  set $p1 [ expr $z_left ]
  set $p2 [ expr $z_right ]

  set z [ expr - $length ]
  while { $z < $length } {
    if { $z < $z_left } {
      set radius [ expr $c1_r - sqrt(  $smoothing_radius*$smoothing_radius - ($z-$c1_z)*($z-$c1_z) ) ] 
      set delta_b [ expr 2*$pi*$res/2/$pi/$smoothing_radius ]
      set sinb [ expr ($z - $c1_z)/$smoothing_radius ]
      set sinbnew [ expr $sinb*cos($delta_b) + sqrt(1-$sinb*$sinb)*sin($delta_b) ]
      set incr_z [ expr $c1_z + $smoothing_radius * $sinbnew - $z ]
      set slope_norm [ expr tan(asin($sinb)) ]
    } elseif  { $z > $z_right } {
      set radius [ expr $c2_r - sqrt(  $smoothing_radius*$smoothing_radius - ($z-$c2_z)*($z-$c2_z) ) ] 
      set delta_b [ expr 2*$pi*$res/2/$pi/$smoothing_radius ]
      set sinb [ expr ($z - $c2_z)/$smoothing_radius ]
      set sinbnew [ expr $sinb*cos($delta_b) + sqrt(1-$sinb*$sinb)*sin($delta_b) ]
      set incr_z [ expr $c2_z + $smoothing_radius * $sinbnew - $z ]
      set slope_norm [ expr tan(asin($sinb)) ]
      if { $incr_z <= 0 } { 
        set z 1000
        continue
      }
    } else {
      set radius [ expr ($z-$c1_z)*$slope + $c1_r -$smoothing_radius/$cosa]  
      set incr_z [expr $res / sqrt(1+$slope*$slope) ]
      set slope_norm $slope
    }
    set n_circle [ expr round( 2*$pi*$radius/$res) ]
    set incr_phi [ expr 2*$pi / $n_circle ]
    set phi 0
    for { set j 0 } { $j < $n_circle } { incr j} {
      set px [ expr $circle_center_x + $radius*cos($phi)*$e_1_x + $radius*sin($phi)*$e_2_x ]
      set py [ expr $circle_center_y + $radius*cos($phi)*$e_1_y + $radius*sin($phi)*$e_2_y ]
      set pz [ expr $circle_center_z + $radius*cos($phi)*$e_1_z + $radius*sin($phi)*$e_2_z ]
      part [ expr $n_induced_charges ] pos $px $py $pz type $type q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] fix 1 1 1
      set nx [ expr -1./sqrt(1+$slope_norm*$slope_norm)*(cos($phi)*$e_1_x + sin($phi)*$e_2_x)+$slope_norm/sqrt(1+$slope_norm*$slope_norm)*$axis_x ]
      set ny [ expr -1./sqrt(1+$slope_norm*$slope_norm)*(cos($phi)*$e_1_y + sin($phi)*$e_2_y)+$slope_norm/sqrt(1+$slope_norm*$slope_norm)*$axis_y ]
      set nz [ expr -1./sqrt(1+$slope_norm*$slope_norm)*(cos($phi)*$e_1_z + sin($phi)*$e_2_z)+$slope_norm/sqrt(1+$slope_norm*$slope_norm)*$axis_z ]
      lappend icc_normals [ list $nx $ny $nz ]
      incr particle_counter
      incr n_induced_charges
      set phi $phi+$incr_phi
    }
    set incr_length_x [ expr $axis_x * $incr_z]
    set incr_length_y [ expr $axis_y * $incr_z]
    set incr_length_z [ expr $axis_z * $incr_z]
    set circle_center_x [ expr $circle_center_x + $incr_length_x ]
    set circle_center_y [ expr $circle_center_y + $incr_length_y ]
    set circle_center_z [ expr $circle_center_z + $incr_length_z ]
    set z [ expr $z + $incr_z ]
  }
  for { set radius [ expr $c1_r + $res ] } { $radius < $max_radius } { set radius [ expr $radius + $res ] } {
    set circle_center_x [ expr $center_x - $axis_x*$length ]
    set circle_center_y [ expr $center_y - $axis_y*$length ]
    set circle_center_z [ expr $center_z - $axis_z*$length ]
    set n_circle [ expr round( 2*$pi*$radius/$res) ]
    set incr_phi [ expr 2*$pi / $n_circle ]
    set phi 0
    for { set j 0 } { $j < $n_circle } { incr j} {
      set px [ expr $circle_center_x + $radius*cos($phi)*$e_1_x + $radius*sin($phi)*$e_2_x ]
      set py [ expr $circle_center_y + $radius*cos($phi)*$e_1_y + $radius*sin($phi)*$e_2_y ]
      set pz [ expr $circle_center_z + $radius*cos($phi)*$e_1_z + $radius*sin($phi)*$e_2_z ]
      if { $px > 0 && $px < $box_l_x && $py > 0 && $py < $box_l_y &&$pz > 0 && $pz < $box_l_z } {  
        part [ expr $n_induced_charges ] pos $px $py $pz type $type q [ expr $sigma*$res*$res +0.1*([ t_random ]-0.5) ] fix 1 1 1 
        set nx [ expr -$axis_x ]
        set ny [ expr -$axis_y ]
        set nz [ expr -$axis_z ]
        lappend icc_normals [ list $nx $ny $nz ]
        incr particle_counter
        incr n_induced_charges
      } 
      set phi $phi+$incr_phi
    }
  }
  for { set radius [ expr $c2_r + $res ] } { $radius < $max_radius } { set radius [ expr $radius + $res ] } {
    set circle_center_x [ expr $center_x + $axis_x*$length ]
    set circle_center_y [ expr $center_y + $axis_y*$length ]
    set circle_center_z [ expr $center_z + $axis_z*$length ]
    set n_circle [ expr round( 2*$pi*$radius/$res) ]
    set incr_phi [ expr 2*$pi / $n_circle ]
    set phi 0
    for { set j 0 } { $j < $n_circle } { incr j} {
      set px [ expr $circle_center_x + $radius*cos($phi)*$e_1_x + $radius*sin($phi)*$e_2_x ]
      set py [ expr $circle_center_y + $radius*cos($phi)*$e_1_y + $radius*sin($phi)*$e_2_y ]
      set pz [ expr $circle_center_z + $radius*cos($phi)*$e_1_z + $radius*sin($phi)*$e_2_z ]
      if { $px > 0 && $px < $box_l_x && $py > 0 && $py < $box_l_y && $pz > 0 && $pz < $box_l_z } {  
        part [ expr $n_induced_charges ] pos $px $py $pz type $type fix 1 1 1
        set nx [ expr +$axis_x ]
        set ny [ expr +$axis_y ]
        set nz [ expr +$axis_z ]
        lappend icc_normals [ list $nx $ny $nz ]
        incr particle_counter
        incr n_induced_charges
      } 
      set phi [ expr $phi+$incr_phi ]
    }
  }
  for { set i 0 } { $i < $particle_counter } { incr i } {
    lappend icc_areas [ expr $res*$res ]
    lappend icc_epsilons $eps
    lappend icc_sigmas $sigma 
  }

}
