setmd box_l 8 8 8 
setmd time_step 1
setmd skin 0.01
electrokinetics agrid 0.05 lb_density 1.0 viscosity 1.0 friction 1.0 T 1.0 bjerrum_length 1.0
electrokinetics 1 density 1 D 0.01 valency 1
electrokinetics pdb-parse polyat.pdb topolat.top
electrokinetics print boundary vtk boundary_at.vtk
electrokinetics 1 print density vtk density_at.vtk
electrokinetics print potential vtk potential_at.vtk
