setmd box_l 10 10 10

setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 4 5 5 q +1 v 0 0 0
part 1 pos 6 5 5 q -1 v 0 0 0
# in real space:
# inter coulomb 1.0 p3m 3 32 5 0.001
# in k-space:
inter coulomb 1.0 p3m gpu 0.0 32 5 1.0

integrate 0
puts [ part 0 print f ]

exit 0
