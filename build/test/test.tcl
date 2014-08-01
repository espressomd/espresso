setmd box_l 20 20 20
setmd time_step 1
setmd skin 0.1
thermostat off

electrokinetics agrid 0.5 lb_density 33.4 viscosity 1 T 1 bjerrum_length 20 friction 1
electrokinetics 1 density 0.06 D 0.01 valency 1 ext_force 0 0 0

for { set i 0} {$i < 100} {incr i} {
  electrokinetics 1 print density vtk data/dens_$i.vtk
  electrokinetics print potential vtk data/pot_$i.vtk
  integrate 10
}
