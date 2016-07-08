# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field outside the slit

set box_x 6
set box_y 6
set width 50

set padding 6
set box_z [expr $width+2*$padding]

setmd box_l $box_x $box_y $box_z

# Set the electrokinetic parameters

set agrid 1.0
set dt 0.2
set kT 1.0
set bjerrum_length 0.7095
set D 0.006075
set valency 1.0
set viscosity_dynamic 79.53
set density_water 26.15
set sigma -0.05
set ext_force 0.1

# Set the simulation parameters

setmd time_step $dt
setmd skin 0.1
thermostat off
set integration_length 200000

# Set up the (LB) electrokinetics fluid

set viscosity_kinematic [expr $viscosity_dynamic/$density_water]
electrokinetics agrid $agrid lb_density $density_water viscosity $viscosity_kinematic friction 1.0 T $kT bjerrum_length $bjerrum_length

# Set up the charged and neutral species

set density_counterions [expr -2.0*double($sigma)/double($width)]
electrokinetics 1 density $density_counterions D $D valency $valency ext_force $ext_force 0 0

# Set up the charged boundaries

electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding-$agrid] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside
electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding+$width] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside

# Set up the walls confining the fluid

electrokinetics boundary charge_density 0.0 wall normal 0 0 1 d $padding 0 0 direction outside
electrokinetics boundary charge_density 0.0 wall normal 0 0 -1 d -[expr $padding+$width] 0 0 direction outside

# Integrate the system

integrate $integration_length

# Output

set fp [open "eof_electrokinetics.dat" "w"]
puts $fp "#position measured_density measured_velocity measured_pressure_xz"

for {set i 0} {$i < [expr $box_z/$agrid]} {incr i} {
  if {[expr $i*$agrid] >= $padding && [expr $i*$agrid] < [expr $box_z - $padding] } {
    set xvalue [expr $i*$agrid - $padding]
    set position [expr $i*$agrid - $padding - $width/2.0 + $agrid/2.0]

    # density
    set measured_density [electrokinetics 1 node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print density]

    # velocity
    set measured_velocity [lindex [electrokinetics node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print velocity] 0]

    # xz component pressure tensor
    set measured_pressure_xz [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 3]

    puts $fp "$position $measured_density $measured_velocity $measured_pressure_xz"
  }
}

close $fp
