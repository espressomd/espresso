# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field outside the slit

set box_x 6
set box_y 6
set width 50

set padding 6
set box_z [expr $width+2*$padding]

# Set the electrokinetic parameters

set agrid 1.0
set temperature 1.0
set bjerrum_length 0.7095
set valency 1.0
set viscosity_dynamic 79.53
set density_water 26.15
set sigma -0.05
set force 0.1

set viscosity_kinematic [expr $viscosity_dynamic/$density_water]
set density_counterions [expr -2.0*double($sigma)/double($width)]

# Calculate the inverse length xi, which is a combination of various
# constants (xi = zeC/2kBT), with C a constant that needs to be
# solved for, or equivalently, xi needs to be solved for

# root finding function

proc solve {xi d bjerrum_length sigma valency } {
  set pi [expr {acos(-1.0)}]
  set el_char 1.0
  return [expr $xi*tan($xi*$d/2.0) + 2.0*$pi*$bjerrum_length*$sigma/($valency*$el_char)]
}

# initial parameters for bisection scheme

set pi [expr {acos(-1.0)}]
set size [expr $pi/(2.0*$width)]

set pnt0 0.0
set pntm [expr $pnt0 + $size]
set pnt1 [expr $pnt0 + 1.9*$size]

# the bisection scheme

set tol 1.0e-08
while { $size > $tol } {

  set val0 [solve $pnt0 $width $bjerrum_length $sigma $valency]
  set val1 [solve $pnt1 $width $bjerrum_length $sigma $valency]
  set valm [solve $pntm $width $bjerrum_length $sigma $valency]

  if { $val0 < 0.0 && $val1 > 0.0 } {
    if { $valm < 0.0 } {
      set pnt0 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt0 + $size]
    } else {
      set pnt1 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt1 - $size]
    }
  } elseif { $val0 > 0.0 && $val1 < 0.0 } {
    if { $valm < 0.0 } {
      set pnt1 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt1 - $size]
    } else {
      set pnt0 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt0 + $size]
    }
  } else {
    error_exit "Bisection method fails:\nTuning of domain boundaries may be required."
  }
}

# obtain the desired xi value

set xi $pntm

# function to calculate the density

proc density {x xi bjerrum_length} {
  set pi [expr {acos(-1.0)}]
  set kb 1.0
  return [expr ($xi*$xi)/(2.0*$pi*$bjerrum_length*cos($xi*$x)*cos($xi*$x))]
}

# function to calculate the velocity

proc velocity {x xi d bjerrum_length force viscosity_kinematic density_water} {
  set pi [expr {acos(-1.0)}]
  return [expr ($force)*log(cos($xi*$x)/cos($xi*$d/2.0))/(2.0*$pi*$bjerrum_length*$viscosity_kinematic*$density_water)]
}

# function to calculate the nonzero component of the pressure tensor

proc pressure_tensor_offdiagonal {x xi bjerrum_length force} {
  set pi [expr {acos(-1.0)}]
  return [expr $force*$xi*tan($xi*$x)/(2.0*$pi*$bjerrum_length)]
}

# function to calculate the hydrostatic pressure

# Technically, the LB simulates a compressible fluid, whiches pressure
# tensor contains an additional term on the diagonal, proportional to 
# the divergence of the velocity. We neglect this (small) contribution.
# The LB pressure tensor also contains the convective acceleration, which
# we neglect here.

proc hydrostatic_pressure {x xi bjerrum_length tensor_entry} {
  return 0.0
}

# Output

set fp [open "eof_analytical.dat" "w"]
puts $fp "#position calculated_density calculated_velocity calculated_pressure_xz"

for {set i 0} {$i < [expr $box_z/$agrid]} {incr i} {

  if {[expr $i*$agrid] >= $padding && [expr $i*$agrid] < [expr $box_z - $padding] } {
    set xvalue [expr $i*$agrid - $padding]
    set position [expr $i*$agrid - $padding - $width/2.0 + $agrid/2.0]

    # density
    set calculated_density [density $position $xi $bjerrum_length]

    # velocity
    set calculated_velocity [velocity $position $xi $width $bjerrum_length $force $viscosity_kinematic $density_water]

    # xz component pressure tensor
    set calculated_pressure_xz [pressure_tensor_offdiagonal $position $xi $bjerrum_length $force]

    puts $fp "$position $calculated_density $calculated_velocity $calculated_pressure_xz"
  }
}

close $fp
