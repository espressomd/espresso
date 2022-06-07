import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import espressomd.shapes


@utx.skipIfMissingFeatures(["LB_WALBERLA"])
class EKNoFlux(ut.TestCase):
    BOX_L = 15.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 50
    RADIUS = 5.

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def test_noflux(self):
        for single_precision in (False, True):
            with self.subTest(single_precision=single_precision):
                self.detail_test_noflux(single_precision=single_precision)

    def detail_test_noflux(self, single_precision: bool):
        """
        Testing the EK noflux boundaries to not leak density outside of a sphere.
        """

        decimal_precision: int = 7 if single_precision else 10

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                   density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0,
                                                   advection=False, friction_coupling=False, ext_efield=[0, 0, 0], single_precision=single_precision)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        center = np.asarray(self.system.box_l / 2, dtype=np.int)

        ekspecies[center[0], center[1], center[2]].density = self.DENSITY

        sphere = espressomd.shapes.Sphere(
            center=self.system.box_l / 2,
            radius=self.RADIUS,
            direction=-1)
        ekspecies.add_boundary_from_shape(
            sphere, [0, 0, 0], espressomd.EKSpecies.FluxBoundary)

        positions = np.empty((*self.system.box_l.astype(np.int), 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, self.system.box_l))
        positions += 0.5

        self.system.integrator.run(self.TIME)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved globally
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, decimal_precision)

        domain_density = simulated_density[np.logical_not(
            ekspecies[:, :, :].is_boundary)]
        # check that the density is kept constant inside the sphere
        np.testing.assert_almost_equal(
            np.sum(domain_density), self.DENSITY, decimal_precision)
        np.testing.assert_array_less(
            0., domain_density, "EK density array contains negative densities!")


if __name__ == "__main__":
    ut.main()
