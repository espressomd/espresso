import espressomd

system = espressomd.System(box_l=3 * [1])

print(system.get_metadata("lennard-jones.pt"))
