from __future__ import print_function
import unittest as ut
import espressomd


@ut.skipIf(not espressomd.has_features(["P3M", "EXTERNAL_FORCES"]),
           "Features not available, skipping test!")
class test_icc(ut.TestCase):

    def runTest(self):
        from espressomd.electrostatics import P3M
        from espressomd.electrostatic_extensions import ICC

        S = espressomd.System(box_l=[1.0, 1.0, 1.0])
        S.seed = S.cell_system.get_state()['n_nodes'] * [1234]
        # Parameters
        box_l = 20.0
        nicc = 10
        q_test = 10.0
        q_dist = 5.0

        # System
        S.box_l = [box_l, box_l, box_l + 5.0]
        S.cell_system.skin = 0.4
        S.time_step = 0.01

        # ICC particles
        nicc_per_electrode = nicc * nicc
        nicc_tot = 2 * nicc_per_electrode
        iccArea = box_l * box_l / nicc_per_electrode

        iccNormals = []
        iccAreas = []
        iccSigmas = []
        iccEpsilons = []

        l = box_l / nicc
        for xi in range(nicc):
            for yi in range(nicc):
                S.part.add(pos=[l * xi, l * yi, 0], q=-0.0001, fix=[1, 1, 1])
                iccNormals.append([0, 0, 1])

        for xi in range(nicc):
            for yi in range(nicc):
                S.part.add(pos=[l * xi, l * yi, box_l],
                           q=0.0001, fix=[1, 1, 1])
                iccNormals.append([0, 0, -1])

        iccAreas.extend([iccArea] * nicc_tot)
        iccSigmas.extend([0] * nicc_tot)
        iccEpsilons.extend([10000000] * nicc_tot)

        # Test Dipole
        b2 = box_l * 0.5
        S.part.add(pos=[b2, b2, b2 - q_dist / 2], q=q_test, fix=[1, 1, 1])
        S.part.add(pos=[b2, b2, b2 + q_dist / 2], q=-q_test, fix=[1, 1, 1])

        # Actors
        p3m = P3M(prefactor=1, mesh=32, cao=7, accuracy=1e-5)
        icc = ICC(
            n_icc=nicc_tot,
            convergence=1e-6,
            relaxation=0.75,
            ext_field=[
                0,
                0,
                0],
            max_iterations=100,
            first_id=0,
            eps_out=1,
            normals=iccNormals,
            areas=iccAreas,
            sigmas=iccSigmas,
            epsilons=iccEpsilons)
        
        S.actors.add(p3m)
        S.actors.add(icc)

        # Run
        S.integrator.run(0)

        # Analyze
        QL = sum(S.part[:nicc_per_electrode].q)
        QR = sum(S.part[nicc_per_electrode:nicc_tot].q)

        testcharge_dipole = q_test * q_dist
        induced_dipole = 0.5 * (abs(QL) + abs(QR)) * box_l

        # Result
        self.assertAlmostEqual(1, induced_dipole / testcharge_dipole, places=4)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
