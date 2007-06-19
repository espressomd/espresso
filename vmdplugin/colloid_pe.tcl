puts "[code_info]"

set vtf_filename "colloid_pe.vtf"

# parameters
set id_colloid 0
set r_colloid 2.0
set n_colloid 5

set id_charged_monomer 1
set id_neutral_monomer 2
set r_monomer 0.5
set n_polymer 5
set l_polymer 30
set q_monomer +1

set box_length 15.0

set n_cycle 50
set n_steps 1000

set time_step 0.001
set skin 0.4
set temp 1.0
set gamma 1.0

set pressure 0.02
set piston_mass 0.0001
set gamma_0 0.5
set gamma_v 0.0001

set fene_epsilon 7.0
set fene_r 2.0
set lj_epsilon 1.0
set lj_sigma 1.0

# derived parameters
set z_colloid [expr ( - $q_monomer * $l_polymer * 0.5) / $n_colloid]
set lj_cutoff [expr $lj_sigma*pow(2.,1./6.)]
set lj_shift [expr 1./4.*$lj_epsilon]

# init simulation
setmd box_l $box_length $box_length $box_length
setmd time_step $time_step
setmd skin $skin
setmd max_num_cells 1000
thermostat langevin $gamma $temp
if {[has_feature NPT]} then {
    integrate set npt_isotropic $pressure $piston_mass 
    thermostat set npt_isotropic $temp $gamma_0 $gamma_v
} else {
    integrate set nvt
}

puts "Simulate PE close to a charged colloid"

puts "Setting global parameters..."
puts "  box_length=[setmd box_l]"
puts "  time_step=[setmd time_step]"
puts "  skin=[setmd skin]"
puts "  integrate=[integrate]"
puts "  thermostat=[thermostat]"

puts "Setting up particles and interactions..."
# steric LJ potentials
inter $id_colloid $id_colloid \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr $r_colloid+$r_colloid - $lj_sigma ]

inter $id_colloid $id_neutral_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr $r_monomer+$r_colloid - $lj_sigma ]

inter $id_neutral_monomer $id_neutral_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr 2*$r_monomer - $lj_sigma ]

inter $id_charged_monomer $id_neutral_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr 2*$r_monomer - $lj_sigma ]

inter $id_charged_monomer $id_charged_monomer \
    lennard-jones \
    $lj_epsilon $lj_sigma \
    $lj_cutoff $lj_shift \
    [ expr 2*$r_monomer - $lj_sigma ]

# bond potential
inter 0 fene $fene_epsilon $fene_r

# create colloids
counterions $n_colloid \
    start 0 charge $z_colloid type $id_colloid

# create the polyelectrolytes
polymer $n_polymer $l_polymer $lj_sigma \
    start [setmd n_part] \
    charge $q_monomer \
    distance 2 \
    types $id_neutral_monomer $id_charged_monomer \
    bond 0 


puts "Writing configurations to $vtf_filename..."
set vtf_file [open $vtf_filename w]
writevsf $vtf_file
writevcf $vtf_file

puts "Warming up..."
set min 0
set cap 10
set warmup_finished 0
while { ! $warmup_finished } {
    # set ljforcecap
    inter ljforcecap $cap
    # integrate a number of steps, e.g. 20
    integrate 20
    for { set pid 0 } { $pid < [setmd n_part] } { incr pid } {
	part $pid v 0 0 0
    }

    set monomer_list [list $id_neutral_monomer $id_charged_monomer ]
    # check the status of the sytem
    set f_cc [expr [analyze mindist $id_colloid $id_colloid] / (2*$r_colloid)]
    set f_mm [expr [analyze mindist $monomer_list $monomer_list] / (2.0*$r_monomer)]
    set f_cm [expr [analyze mindist $id_colloid $monomer_list] / (0.8*($r_colloid+$r_monomer))]
    puts "f=($f_cc, $f_mm, $f_cm)"
    set warmup_finished [expr $f_cc > 0.8 && $f_mm > 0.8 && $f_cm > 0.8]

    # this is a shortcut for 'set cap [expr $cap+10]'
    incr cap 10
}
puts "Warmup finished."
writevcf $vtf_file

# turn off the ljforcecap, which is done by setting the 
# force cap value to zero:
inter ljforcecap 0
# reset velocities
for { set pid 0 } { $pid < [setmd n_part] } { incr pid } {
    part $pid v 0 0 0
}
# reset simulation time
setmd time 0.0

puts "Tuning P3M..."
#puts [inter coulomb 1 p3m tunev2 accuracy 0.01]
#inter coulomb 1 p3m 13.9629 8 3 0.127621 0.000832468

puts "Starting simulation..."
set i 0 
while { $i<$n_cycle } {
    puts "$i [setmd box_l]"
    integrate $n_steps

    writevcf $vtf_file

    incr i
}
close $vtf_file

puts "Simulation finished."

