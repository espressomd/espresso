
setmd box_l 10 10 10


#dielectric_sphere center 5 5 5 radius 3 res 1 eps 0.2 
#global n_induced_charges icc_areas icc_normals icc_epsilons icc_sigmas
#puts [ llength $icc_areas ] 
#puts [ llength $icc_normals ]
#puts [ llength $icc_epsilons ] 
#puts [ llength $icc_sigmas ]
#dielectric_sphere center 10 10 10 radius 3 res 1 eps 0.2 
#
#puts [ llength $icc_areas ] 
#puts [ llength $icc_normals ]
#puts [ llength $icc_epsilons ] 
#puts [ llength $icc_sigmas ]
#puts "npart [ setmd n_part ]"
#
#dielectric_cylinder center 13 14 15 axis 3 3 3 radius 4 length 10 res 1 eps 0.2 covers 1
#puts "areas [ llength $icc_areas ] "
#puts "norms [ llength $icc_normals ]"
#puts "epsil [ llength $icc_epsilons ]" 
#puts "sigma [ llength $icc_sigmas ]"
#puts [ setmd n_part ]
#
#dielectric_wall dist 1. normal 1 0 0 res 1
#dielectric_wall dist 1. normal 0 1 0 res 1
#dielectric_wall dist 1. normal 0 0 1 res 1
dielectric_pore  center 10 10 10 axis 1 0 0 radii 4.5 8.5 length 11. res 0.1 smoothing_radius 1.5

set ofile [ open "diel.dat" "w" ]
for { set i 0 } { $i < $n_induced_charges } { incr i } {
  puts $ofile "[ part $i print pos ] [ lindex $icc_normals $i ]" 
}
#puts "n_induced_charges $n_induced_charges"
#puts "areas [ llength $icc_areas ] "
#puts "norms [ llength $icc_normals ]"
#puts "epsil [ llength $icc_epsilons ]" 
#puts "sigma [ llength $icc_sigmas ]"
#puts [ setmd n_part ]

constraint pore center 10 10 10 axis 1 0 0 radii 5. 9. length 5  type 0

set distfile [open "dist.dat" "w" ]
for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
  set pos [ part $i print pos ]
  puts $distfile "[ lindex $pos 0 ] [ constraint mindist_position [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ] ]"
}

#prepare_vmd_connection "test" 3000
