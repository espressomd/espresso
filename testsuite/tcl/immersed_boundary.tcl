source "tests_common.tcl"

require_feature "LB"
require_feature "LB_BOUNDARIES"
require_feature "IMMERSED_BOUNDARY"

puts "----------------------------------------"
puts "- Testcase immersed_boundary.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

if { [catch {
  ###### Setup #########
  set boxX 12
  set boxY 12
  set boxZ 20
  setmd box_l $boxX $boxY $boxZ

  set gridDist 1
  set nu 1
  set rho 1
  set mu 1
  set f 1e-3
  
  # time step: use the standard LBM time step
  set taurelax 1
  set dt [expr ($taurelax-0.5)*$gridDist*$gridDist/(3*$nu)]
  
  setmd time_step $dt
  setmd skin 0.1
  lbfluid agrid $gridDist dens $rho visc $nu tau $dt ext_force $f 0 0
  
  # set up boundaries
  lbboundary wall dist 0.5               normal 0 0  1 velocity 0 0 0
  lbboundary wall dist [expr -$boxZ + 0.5] normal 0 0 -1 velocity 0 0 0
  
  # place the virtual particle
  part 0 pos 0 0 [expr $boxZ/2] type 0 virtual 1
  
  ####### integration ######
  
  set numSteps 10000
  integrate $numSteps
  
  ###### check final position, should flow along with the center velocity ####
  
  ### this is what would be expected, however: (1) the profile is not perfectly Poiseuille and (2) it takes some time for the flow to establish
  #set umax [expr ($boxZ-1)*($boxZ-1)/(8*$mu)*$f - 0.005]
  #set xExpect [expr $umax * $dt * $numSteps]
  # therefore we set the value by hand
  set xExpect 66.00
  
  set pos [part 0 print pos]
  set xSim [lindex $pos 0]
  
  set diff [expr abs($xExpect-$xSim)]
  
  if { $diff > 0.01 } {
    error "Immersed boundary virtual particle is not correctly advected by the fluid."
  }
} res ] } {
    error_exit $res
}  
exit 0
  
