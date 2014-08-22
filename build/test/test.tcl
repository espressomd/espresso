setmd box_l 20 20 20
setmd time_step 1
setmd skin 0.1
thermostat off

#electrokinetics agrid 0.5 lb_density 33.4 viscosity 1 T 1 bjerrum_length 0.7 friction 1
#electrokinetics 1 density 0.06 D 0.01 valency 0 ext_force 1 0 0
electrokinetics agrid 1 lb_density 10 viscosity 1 T 1 bjerrum_length 0.7 friction 1 use_nonlinear_stencil 0
electrokinetics 1 density 0.01 D 0.01 valency 0 ext_force -0.577 0.577 0.577

#integrate 10000
#exit

for { set i 0} {$i < 100} {incr i} {
  puts "step $i"
  electrokinetics 1 print density vtk data/dens_$i.vtk
  electrokinetics print potential vtk data/pot_$i.vtk
  electrokinetics print velocity vtk data/vel_$i.vtk
  electrokinetics 1 print flux vtk data/flux_$i.vtk
  electrokinetics print lbforce vtk data/force_$i.vtk
  #lbfluid print vtk force data/force_$i.vtk

  integrate 1
}
