
source "tests_common.tcl"

require_feature "CONSTRAINTS"

puts "---------------------------------------------------------------"
puts "- Testcase dielectric_test.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 20 20 20

set dist 0.5

puts "Creating Wall 1"
dielectric_wall normal 1. 0 0 dist 0.5 res 1.2
constraint wall normal 1. 0 0 dist -0.0 type 0
puts "Done. Now we have $n_induced_charges induced charges"
puts "Creating Wall 2"
dielectric_wall normal -1. 0 0 dist -19 res 1.2
constraint wall normal -1. 0 0 dist -19.5 type 0
puts "Done. Now we have $n_induced_charges induced charges"


puts "Creating sphere"
set r 3
dielectric_sphere center 14 10 10 radius $r res 0.5 eps 0.2 
constraint sphere center 14 10 10 radius [ expr $r-$dist ] direction 1 type 0
puts "Done. Now we have $n_induced_charges induced charges"

puts "Creating pore"
set r1 5. 
set r2 8.
set rs 2.5
set l 4.
dielectric_pore  center 10 10 10 axis 1 0 0 radii $r1 $r2 length $l res 1.9 smoothing_radius $rs eps 0.2

constraint pore center 10 10 10 axis 1 0 0 radii [ expr $r1 + $dist ] [ expr $r2 + $dist ] length [ expr $l - $dist ] type 0 smoothing_radius [ expr $rs - $dist ]
puts "Done. Now we have $n_induced_charges induced charges"


puts "The cylinder test is still disable because it is not fully consistent"
puts "with the cylinder constraint"
#puts "Creating cylinder"
#set r 2
#set l 3
#dielectric_cylinder center 10 10 10   axis 1 0 0 radius $r length $l res .25 covers 1
#constraint cylinder center 10 10 10   axis 1 0 0 radius [ expr $r - $dist ] length [ expr $l - $dist ] type 0 direction 1
#puts "Done. Now we have $n_induced_charges induced charges"

puts "Checking all distances and normal vectors"
for { set i 0 } { $i < $n_induced_charges } { incr i } {
  set pos [ part $i print pos ]
  set dv [ constraint mindist_position_vec [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
  set n  [ lindex $icc_normals $i ] 
  set diff [ expr [ lindex $dv 0 ] /$dist - [ lindex $n 0 ]] 
  if { $diff > 1e-6 } {
    puts "Everything is wrong at id $i!"
    puts "pos is $pos dist_vec is $dv normal is $n"
    exit 1
  }
}
puts "Test is fully passed."

puts "n_induced_charges $n_induced_charges"
puts "areas [ llength $icc_areas ] "
puts "norms [ llength $icc_normals ]"
puts "epsil [ llength $icc_epsilons ]" 
puts "sigma [ llength $icc_sigmas ]"
puts [ setmd n_part ]
exit 0

##### We can enable this to make visual checks of everything!

set distfile [ open "dist.dat" "w" ]
for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
  set pos [ part $i print pos ]
  set dv [ constraint mindist_position_vec [ lindex $pos 0 ] [ lindex $pos 1 ] [ lindex $pos 2 ]   ]
  puts $distfile "$pos $dv" 
}
close $distfile

#prepare_vmd_connection "test" 3000
#after 200000

