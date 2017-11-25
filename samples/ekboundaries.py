from espressomd import System, shapes, electrokinetics
import sys
system = System()
system.box_l = [10, 10, 10]
system.cell_system.skin = 0.4
system.time_step = 0.1


ek = electrokinetics.Electrokinetics(
    lb_density=1, friction=1, agrid=1, viscosity=1, T=1, prefactor=1)

pos = electrokinetics.Species(
    density=0.05, D=0.1, valency=1, ext_force=[0, 0, 1.])
neg = electrokinetics.Species(
    density=0.05, D=0.1, valency=-1, ext_force=[0, 0, -1.])
ek.add_species(pos)
ek.add_species(neg)
system.actors.add(ek)

print(ek.get_params())
print(pos.get_params())
print(neg.get_params())
print(pos[5, 5, 5].density)


ek_wall_left = electrokinetics.EKBoundary(
    shape=shapes.Wall(dist=1, normal=[1, 0, 0]), charge_density=-0.01)
ek_wall_right = electrokinetics.EKBoundary(
    shape=shapes.Wall(dist=-9, normal=[-1, 0, 0]), charge_density=0.01)
system.ekboundaries.add(ek_wall_left)
system.ekboundaries.add(ek_wall_right)

for i in range(1000):
    system.integrator.run(100)
    sys.stdout.write("\rIntegrating: %03i" % i)
    sys.stdout.flush()

    pos.print_vtk_density("ek/pos_dens_%i.vtk" % i)
    neg.print_vtk_density("ek/neg_dens_%i.vtk" % i)
    pos.print_vtk_flux("ek/pos_flux_%i.vtk" % i)
    neg.print_vtk_flux("ek/neg_flux_%i.vtk" % i)
    ek.print_vtk_velocity("ek/ekv_%i.vtk" % i)
    ek.print_vtk_boundary("ek/ekb_%i.vtk" % i)
