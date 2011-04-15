setmd box_l 10 10 10

setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m 3 32 5 0.001
# in k-space:
inter coulomb 1.0 p3m 0.0 32 5 1.0

integrate 0
puts [ part 0 print f ]
set f1 [ lindex [part 0 print f] 0 ]

setmd box_l 20 10 10

invalidate_system 

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
part 2 pos 14 5 5 q +1 v 0 0 0
part 3 pos 16 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m 3 32 5 0.001
# in k-space:
inter coulomb 1.0 p3m 0.0 64 5 1.0
integrate 0

puts [ part 0 print f ]
set f2 [ lindex [part 0 print f] 0 ]

if { abs( $f1 - $f2 )/$f1 < 1e-5 } {
  puts OK
} else {
  error "P3M noncubic test failed"
}

integrate 0
