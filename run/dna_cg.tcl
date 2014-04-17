setmd box_l 80 80 80
setmd time_step 0.2
setmd skin 0.1
thermostat off

electrokinetics agrid 1 lb_density 33.0 viscosity 1.0 friction 1.0 T 1.0 bjerrum_length 0.7
electrokinetics 1 density 0.001 D 0.3 valency 1 ext_force 0 0 0.001
electrokinetics pdb-parse polycg.pdb topolcg.top

electrokinetics print boundary vtk data/boundary_cg.vtk

for {set i 0} {1} {incr i} {
  puts $i
  integrate 10
  electrokinetics 1 print density vtk data/density_cg_$i.vtk
  electrokinetics 1 print flux vtk data/flux_cg_$i.vtk
  electrokinetics print velocity vtk data/velocity_cg_$i.vtk
  electrokinetics print potential vtk data/potential_cg_$i.vtk
}
