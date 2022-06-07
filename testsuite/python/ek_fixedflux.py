import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.shapes


@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class EKFixedFlux(ut.TestCase):
    BOX_L = 5.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 10

    INFLOW_FLUX = 0.1

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def test_inflow(self):
        for single_precision in (False, True):
            with self.subTest(single_precision=single_precision):
                self.detail_test_inflow(single_precision=single_precision)

    def detail_test_inflow(self, single_precision: bool):
        """
        Testing the EK fixed flux boundaries to test the fixed inflow into a non-periodic box.
        """

        decimal_precision: int = 5 if single_precision else 7

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                   density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
                                                   valency=0.0, advection=False, friction_coupling=False,
                                                   ext_efield=[0, 0, 0], single_precision=single_precision)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        ekspecies[1:-1, 1:-1, 1:-1].density = self.DENSITY

        ekspecies[:, :, 0].flux_boundary = espressomd.EKSpecies.FluxBoundary([
            0, 0, 0])
        ekspecies[:, :, -
                  1].flux_boundary = espressomd.EKSpecies.FluxBoundary([0, 0, 0])
        ekspecies[:, 0, :].flux_boundary = espressomd.EKSpecies.FluxBoundary([
            0, 0, 0])
        ekspecies[:, -1,
                  :].flux_boundary = espressomd.EKSpecies.FluxBoundary([0, 0, 0])
        ekspecies[0, :, :].flux_boundary = espressomd.EKSpecies.FluxBoundary([
            0, 0, 0])
        ekspecies[-1, :,
                  :].flux_boundary = espressomd.EKSpecies.FluxBoundary([0, 0, 0])

        # set fixed flux in +z-direction
        ekspecies[:, :, 4].flux_boundary = espressomd.EKSpecies.FluxBoundary(
            [0, 0, -self.INFLOW_FLUX])
        additional_center_flux = 3 * self.INFLOW_FLUX
        ekspecies[int(self.BOX_L //
                      2), int(self.BOX_L //
                              2), 4].flux_boundary = espressomd.EKSpecies.FluxBoundary([0, 0, -
                                                                                        self.INFLOW_FLUX -
                                                                                        additional_center_flux])

        # check density before integration
        expected_initial_density = self.DENSITY * (self.BOX_L - 2)**3

        np.testing.assert_almost_equal(actual=np.sum(
            ekspecies[1:-1, 1:-1, 1:-1].density), desired=expected_initial_density, decimal=decimal_precision)

        self.system.integrator.run(self.TIME)

        # check that density has pushed into domain
        inflow_area = (ekspecies.shape[0] - 2) * (ekspecies.shape[1] - 2)
        expected_end_density = expected_initial_density + \
            (self.INFLOW_FLUX * inflow_area + additional_center_flux) * self.TIME

        np.testing.assert_almost_equal(actual=np.sum(
            ekspecies[1:-1, 1:-1, 1:-1].density), desired=expected_end_density, decimal=decimal_precision)


if __name__ == "__main__":
    ut.main()
