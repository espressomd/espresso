#
# Copyright (C) 2022-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy as np
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.electrokinetics


def poiseuille_flow(z, H, ext_force_density, dyn_visc):
    """
    Analytical solution for planar Poiseuille flow.

    Parameters
    ----------
    z : :obj:`float`
        Distance to the mid plane of the channel.
    H : :obj:`float`
        Distance between the boundaries.
    ext_force_density : :obj:`float`
        Force density on the fluid normal to the boundaries.
    dyn_visc : :obj:`float`
        Dynamic viscosity of the LB fluid.

    """
    return ext_force_density * 1. / (2 * dyn_visc) * (H**2.0 / 4.0 - z**2.0)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class FrictionCoupling:
    BOX_L = [28.65, 8.595, 8.595]
    AGRID = 2.865
    TIMESTEPS = 2000
    TAU = 4.483

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.lb = None
        self.system.ekcontainer = None

    def test_friction_coupling(self):
        """
        Use the species as an homogeneous external force on the fluid to check
        the friction coupling against the analytical solution for the planar
        Poiseuille flow.
        """
        kT = 9.8765
        diffusion_coefficient = 0.68
        fluid_density = 1.9
        valency = 0.8
        external_electric_field = np.asarray([0.0, 0.0, 0.01])

        visc = 1. / 6.
        density = 0.0006

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        num_species = 2
        ekspecies = []
        for _ in range(num_species):
            ekspecies.append(espressomd.electrokinetics.EKSpecies(
                lattice=lattice, density=density / num_species, kT=kT, valency=valency,
                diffusion=diffusion_coefficient, friction_coupling=True,
                advection=False, ext_efield=external_electric_field,
                tau=self.TAU, **self.lattice_params))

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU, **self.lattice_params)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        for specie in ekspecies:
            self.system.ekcontainer.add(specie)

        lb_fluid = espressomd.lb.LBFluid(
            lattice=lattice, density=fluid_density, kinematic_viscosity=visc,
            tau=self.TAU, **self.lattice_params)
        self.system.lb = lb_fluid

        wall_top = espressomd.shapes.Wall(normal=[1, 0, 0], dist=self.AGRID)
        wall_bottom = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.BOX_L[0] - self.AGRID))
        lb_fluid.add_boundary_from_shape(wall_top, [0, 0, 0])
        lb_fluid.add_boundary_from_shape(wall_bottom, [0, 0, 0])

        self.system.integrator.run(self.TIMESTEPS)

        velocity_simulation = np.copy(lb_fluid[1:-1, 1, 1].velocity)[:, 2]
        z = (np.arange(self.BOX_L[0] / self.AGRID) +
             0.5) * self.AGRID - self.BOX_L[0] / 2
        # remove the boundary layers
        z = z[1:-1]
        velocity_analytic = poiseuille_flow(
            z, self.BOX_L[0] - 2 * self.AGRID, external_electric_field[2] * valency * density, visc * fluid_density)

        atol = 5e-9 if self.lattice_params["single_precision"] else 8e-16
        np.testing.assert_allclose(
            velocity_simulation, velocity_analytic, atol=atol, rtol=0.0)


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKTestWalberlaDoublePrecisionCPU(FrictionCoupling, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKTestWalberlaSinglePrecisionCPU(FrictionCoupling, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKTestWalberlaDoublePrecisionGPU(FrictionCoupling, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(FrictionCoupling, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
