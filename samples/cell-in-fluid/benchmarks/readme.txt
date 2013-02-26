note1: unless edited otherwise, the script runbenchmarks will run all the
following benchmarks

note2: it is necessary to edit the path in line 2 of runbenchmarks to
point to your working version of espresso

note3: the script runbenchmarks assumes that you are located in the
directory benchmarks  

benchmarks1.tcl works with tetrahedron
it runs the espresso simulation 6 times in a row with the combination
of parameters specified in runbenchmarks. these can be edited, however
note, that after the runs, it compares the outputs with the
precomputed benchmark values stored in folder b_results and notifies
the user of discrepancies.
this script is meant to check the correct behavior of elastic forces
and looks at the position and velocity of all 4 tetrahedron vertices 
at given time intervals 

benchmarks2.tcl works with sphere
(for now) it runs the espresso simulation once 
it simulates pulling of the sphere by force to reach given velocity
then the sphere is released and slows down
the fluid is stationary through the simulation
this script is meant to check the correct interaction between fluid
and object

benchmarks3.tcl works with red blood cell //not prepared yet





 