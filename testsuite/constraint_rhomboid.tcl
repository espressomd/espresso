setmd skin 0.1
setmd box_l 15 15 15
setmd time_step 0.01
setmd periodic 0 0 0
thermostat off

constraint rhomboid corner 5 5 5 a 5 0 0 b 0 5 0 c 0 0 5 direction outside type 0

part 0 pos 7.5 7.5 0 v 0 0 1 type 0
part 1 pos 0 7.5 7.5 v 1 0 0 type 0
part 2 pos 7.5 0 7.5 v 0 1 0 type 0
part 3 pos 7.5 7.5 15 v 0 0 -1 type 0
part 4 pos 15 7.5 7.5 v -1 0 0 type 0
part 5 pos 7.5 15 7.5 v 0 -1 0 type 0

part 6 pos 0 0 0 v 1 1 1 type 0
part 7 pos 0 15 0 v 1 -1 1 type 0
part 8 pos 0 15 15 v 1 -1 -1 type 0
part 9 pos 0 0 15 v 1 1 -1 type 0
part 10 pos 15 0 0 v -1 1 1 type 0
part 11 pos 15 15 0 v -1 -1 1 type 0
part 12 pos 15 15 15 v -1 -1 -1 type 0
part 13 pos 15 0 15 v -1 1 -1 type 0

part 14 pos 7.5 0 0 v 0 1 1 type 0
part 15 pos 0 7.5 0 v 1 0 1 type 0
part 16 pos 0 0 7.5 v 1 1 0 type 0
part 17 pos 7.5 15 15 v 0 -1 -1 type 0
part 18 pos 15 7.5 15 v -1 0 -1 type 0
part 19 pos 15 15 7.5 v -1 -1 0 type 0
part 20 pos 7.5 0 15 v 0 1 -1 type 0
part 21 pos 0 7.5 15 v 1 0 -1 type 0
part 22 pos 15 0 7.5 v -1 1 0 type 0
part 23 pos 0 15 7.5 v 1 -1 0 type 0
part 24 pos 7.5 15 0 v 0 -1 1 type 0
part 25 pos 15 7.5 0 v -1 0 1 type 0

inter 0 0 lennard-jones 1. 1. 1.1225 0.25 0

prepare_vmd_connection "vmd" 2000 1 1

while {1} {
  integrate 1
  imd positions
  after 5
}
