## poisseuille.tcl
## 
## A script to simulate planar Poisseuille flow in Espresso

## set fake MD parameters so that Espresso does 
## not complain
setmd time_step 0.01
setmd skin 0.2

## Set the box size
setmd box_l 16. 16. 16. 

## Set up the LB fluid
#lbfluid agrid 1. dens 10. visc 1. tau 0.1 gamma_odd -0.7 gamma_even -0.7 ext_force 0. 0.1 0.
lbfluid gpu agrid 1. dens 1. visc 1.0 tau 0.01 ext_force 0. 0.001 0.
## set the temperature to a nonzero value to activate fluctuations
thermostat lb 0.

## Set LB boundaries
lbboundary wall dist 1.5 normal 1. 0. 0.
lbboundary wall dist -14.5 normal -1. 0. 0.

## Perform enough iterations until the flow profile
## is stable (1000 LB updates):
integrate 100000

## Part of the solution
set ofile [ open "fluid_velocity.dat" "w" ]
for { set i 0 } { $i < 16 } { incr i } {
	puts $ofile "$i [lindex [lbnode $i 0 0 print u] 1]"
}
close $ofile

