setmd time_step 0.14
setmd skin 0.1
setmd box_l 20 20 20
thermostat off

electrokinetics agrid 1 lb_density 33 viscosity 1 T 1 bjerrum_length 0.7
electrokinetics 1 density 0.1 D 0.3 valency 0
electrokinetics 2 density 0.1 D 0.3 valency 0

for {set i 0} {$i < 100} {incr i} {
    electrokinetics 1 print density vtk data/dens1_$i.vtk
    electrokinetics 2 print density vtk data/dens2_$i.vtk
    integrate 1
}

#electrokinetics pdb-parse cg_dna.pdb cg_dna.itp

#electrokinetics print vtk boundary boundary.vtk
#electrokinetics print vtk potential potential.vtk
