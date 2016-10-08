from espressomd import System, electrokinetics

# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field outside the slit

box_x = 6
box_y = 6
width = 50

padding = 6
box_z = width + 2*padding

system = System()
system.box_l = [box_x, box_y, box_z]

# Set the electrokinetic parameters

agrid = 1.0
dt = 0.2
kT = 1.0
bjerrum_length = 0.7095
D = 0.006075
valency = 1.0
viscosity_dynamic = 79.53
density_water = 26.15
sigma = -0.05
ext_force = 0.1

# Set the simulation parameters

system.time_step = dt
system.cell_system.skin = 0.2
system.thermostat.turn_off()
integration_length = 2e5

# Set up the (LB) electrokinetics fluid
viscosity_kinematic = viscosity_dynamic / density_water
ek = electrokinetics.Electrokinetics(agrid = agrid, lb_density = density_water,
                                     viscosity = viscosity_kinematic, friction = 1.0,
                                     T = kT, bjerrum_length = bjerrum_length)
system.actors.add(ek)

# TODO
# # Set up the charged and neutral species

# set density_counterions [expr -2.0*double($sigma)/double($width)]
# electrokinetics 1 density $density_counterions D $D valency $valency ext_force $ext_force 0 0

# # Set up the charged boundaries

# electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding-$agrid] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside
# electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding+$width] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside

# # Set up the walls confining the fluid

# electrokinetics boundary charge_density 0.0 wall normal 0 0 1 d $padding 0 0 direction outside
# electrokinetics boundary charge_density 0.0 wall normal 0 0 -1 d -[expr $padding+$width] 0 0 direction outside


# Integrate the system
system.integrator.run(integration_length)

# Output
position_list = []
density_list = []
velocity_list = []
pressure_xz_list = []

for i in range(int(box_z/agrid)):
    if (i*agrid >= padding) and (i*agrid < box_z - padding):
        position = i*agrid - padding - width/2.0 + agrid/2.0
        position_list.append(position)

        # TODO
        # # density
        # set measured_density [electrokinetics 1 node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print density]

        # # velocity
        # set measured_velocity [lindex [electrokinetics node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print velocity] 0]

        # # xz component pressure tensor
        # set measured_pressure_xz [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 3]




np.savetxt("eof_analytical.dat", np.column_stack((position_list, density_list, velocity_list, pressure_xz_list)), header="#position calculated_density calculated_velocity calculated_pressure_xz")
