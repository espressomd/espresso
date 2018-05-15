import numpy as np
import matplotlib.pyplot as plt

import espressomd
import espressomd.lb
import espressomd.observables
import espressomd.shapes
import espressomd.lbboundaries
import espressomd.accumulators


system = espressomd.System(box_l=[10.0, 10.0, 5.0])
system.time_step = 0.01
system.cell_system.skin = 0.4

lb_fluid = espressomd.lb.LBFluidGPU(agrid=1.0, fric=1.0, dens=1.0, visc=1.0, tau=0.01, ext_force=[0, 0, 0.15])
system.actors.add(lb_fluid)
system.thermostat.set_lb(kT=1.0)
fluid_obs = espressomd.observables.CylindricalLBVelocityProfile(
        center = [5.0, 5.0, 0.0],
        axis = 'z',
        n_r_bins = 100,
        n_phi_bins = 1,
        n_z_bins = 1,
        min_r = 0.0,
        max_r = 4.0,
        min_phi = -np.pi,
        max_phi = np.pi,
        min_z = 0.0,
        max_z = 10.0,
        sampling_delta_x = 0.05,
        sampling_delta_y = 0.05,
        sampling_delta_z = 1.0)
cylinder_shape = espressomd.shapes.Cylinder(
        center = [5.0, 5.0, 5.0],
        axis = [0, 0, 1],
        direction = -1,
        radius = 4.0,
        length = 20.0)
cylinder_boundary = espressomd.lbboundaries.LBBoundary(shape=cylinder_shape)
system.lbboundaries.add(cylinder_boundary)
system.integrator.run(5000)


accumulator = espressomd.accumulators.MeanVarianceCalculator(obs=fluid_obs)
system.auto_update_accumulators.add(accumulator)
system.integrator.run(5000)

lb_fluid_profile = accumulator.get_mean()
lb_fluid_profile = np.reshape(lb_fluid_profile, (100, 1, 1, 3))

def poiseuille_flow(r, R, ext_force):
    return ext_force * 1./4 * (R**2.0-r**2.0)


# Please note that due to symmetry and interpolation a plateau is seen near r=0.
n_bins = len(lb_fluid_profile[:, 0, 0, 2])
r_max = 4.0
r = np.linspace(0.0, r_max, n_bins)
plt.plot(r, lb_fluid_profile[:, 0, 0, 2], label='LB profile')
plt.plot(r, poiseuille_flow(r, r_max, 0.15), label='analytical solution')
plt.show()
