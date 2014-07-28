set box_x 20 
set box_y 20
set width 20
set padding 6
set box_z [expr $width+2*$padding]
setmd box_l $box_x $box_y $box_z
setmd time_step 0.14 
setmd skin 0.1
thermostat off

electrokinetics agrid 1 lb_density 33 viscosity 1 T 1 bjerrum_length 0.7

electrokinetics 1 density 0.006 D 0.01 valency 1.0
electrokinetics 2 density 0.006 D 0.01 valency -1.0

#electrokinetics boundary charge_density 0.05 rhomboid corner 0 0 [expr $padding-1] c 0 0 1 b $box_x 0 0 a 0 $box_y 0 direction outside
#electrokinetics boundary charge_density -0.05 rhomboid corner 0 0 [expr $padding+$width] c 0 0 1 b $box_x 0 0 a 0 $box_y 0 direction outside

#electrokinetics boundary charge_density 0.0 wall normal 0 0 1 d $padding 0 0 direction outside
#electrokinetics boundary charge_density 0.0 wall normal 0 0 -1 d -[expr $padding+$width] 0 0 direction outside

for {set i 0} {$i < 100} {incr i} {
    electrokinetics 1 print density vtk data/dens1_$i.vtk
    electrokinetics 2 print density vtk data/dens2_$i.vtk
    electrokinetics print potential vtk data/potential_$i.vtk
    integrate 100
}
