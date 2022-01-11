import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.shapes


@utx.skipIfMissingFeatures(["EK_WALBERLA"])
class EKFixedFlux(ut.TestCase):
    BOX_L = 5.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 50

    INFLOW_FLUX = 0.1

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def test_inflow(self):
        """
        Testing the EK fixed flux boundaries to test the fixed inflow into a non-periodic box.
        """

        lattice = espressomd.lb.LatticeWalberla(
            box_size=self.system.box_l, n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                   density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
                                                   valency=0.0, advection=False, friction_coupling=False,
                                                   ext_efield=[0, 0, 0])

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        ekspecies[1:-1, 1:-1, 1:-1].density = self.DENSITY

        offset = 1.0
        # set NoFlux boundary condition on every wall except in +z-direction
        wall_params = [(offset, [1, 0, 0]),
                       (offset, [0, 1, 0]),
                       (offset, [0, 0, 1]),
                       (-(self.BOX_L - offset), [-1, 0, 0]),
                       (-(self.BOX_L - offset), [0, -1, 0])]
        for dist, normal in wall_params:
            wall = espressomd.shapes.Wall(dist=dist, normal=normal)
            ekspecies.add_boundary_from_shape(
                wall, [0, 0, 0], espressomd.EKSpecies.FluxBoundary)

        # set fixed flux in +z-direction
        wall = espressomd.shapes.Wall(
            dist=-(self.BOX_L - offset), normal=[0, 0, -1])

        ekspecies.add_boundary_from_shape(
            wall, [0, 0, self.INFLOW_FLUX], espressomd.EKSpecies.FluxBoundary)

        # check density before integration
        expected_initial_density = self.DENSITY * (self.BOX_L - 2 * offset)**3

        np.testing.assert_almost_equal(actual=np.sum(
            ekspecies[1:-1, 1:-1, 1:-1].density), desired=expected_initial_density)

        self.system.integrator.run(self.TIME)

        # check that density has pushed into domain
        expected_end_density = expected_initial_density + \
            self.INFLOW_FLUX * (self.BOX_L - 2 * offset)**2 * self.TIME

        np.testing.assert_almost_equal(actual=np.sum(
            ekspecies[1:-1, 1:-1, 1:-1].density), desired=expected_end_density)


if __name__ == "__main__":
    ut.main()
