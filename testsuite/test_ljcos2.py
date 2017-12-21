# This is a unit-test for Lennard Jones Cosine2 non-bonded interaction

from __future__ import print_function
import unittest
import espressomd
import numpy as np


@unittest.skipIf(not espressomd.has_features(["LJCOS2"]),
                 "Features not available, skipping test!")
class test_ljCos2(unittest.TestCase):

    test_file = 'data/test_ljcos2.dat'

    # Test Parameters
    epsilon = 1E-4
    numPart = 400

    # System Parameters
    system = espressomd.System()
    system.time_step = 0.001
    system.cell_system.skin = 0.
    system.thermostat.turn_off()

    # Data Containers
    data = []
    prevParticleData = []
    curEnergy = {}
    curPressure = {}
    prevTotalEnergy = 0.
    prevTotalPressure = 0.

    @classmethod
    def read(cls):
        cls.prevTotalEnergy, cls.prevTotalPressure = np.genfromtxt(
            cls.test_file, skip_header=2, max_rows=1, unpack=True)
        cls.prevParticleData = np.genfromtxt(
            cls.test_file, skip_header=8, dtype=[int] * 2 + [float] * 6)
        cls.system.box_l = np.genfromtxt(
            cls.test_file, skip_header=5, max_rows=1, usecols=(1, 2, 3))
        for p in cls.prevParticleData:
            cls.system.part.add(id=p[0], type=p[1], pos=(p[2], p[3], p[4]))

    @classmethod
    def write(cls):
        energy = cls.system.analysis.energy()
        pressure = cls.system.analysis.pressure()
        file_w = open(cls.test_file, 'w')
        file_w.write('\nTotalEnergy\tTotalPressure\n')
        file_w.write('\n%f\t%f\n' % (energy['total'], pressure['total']))
        file_w.write('\nbox_l: %f %f %f\n' %
                     (system.box_l[0], system.box_l[1], system.box_l[2]))
        file_w.write('\nPid\tType\tPosition\tVelocity\tForce\n')

        for i in range(numPart):
            particle = cls.system.part[i]
            file_w.write('\n%d\t%d\t%f %f %f\t%f %f %f' % (particle.id, particle.type,
                                                           particle.pos[0], particle.pos[1], particle.pos[2],
                                                           particle.f[0], particle.f[1], particle.f[2]))

    @classmethod
    def time_translation_dt(cls):
        # Set up LennardJonesCos2
        cls.system.non_bonded_inter[0, 0].lennard_jones_cos2.set_params(
            epsilon=1., sigma=1.,
            width=1.6, offset=0.)
        cls.system.non_bonded_inter[0, 1].lennard_jones_cos2.set_params(
            epsilon=1., sigma=1.,
            width=1.1, offset=0.)
        cls.system.non_bonded_inter[1, 1].lennard_jones_cos2.set_params(
            epsilon=1., sigma=1.,
            width=1.5, offset=0.)
        cls.system.integrator.run(0)

    @classmethod
    def setUpClass(cls):
        cls.read()
        cls.time_translation_dt()
        # After time translation
        cls.curEnergy = cls.system.analysis.energy()
        cls.curPressure = cls.system.analysis.pressure()

    def test_energy(self):
        totalEnergy = test_ljCos2.curEnergy['total']
        nonbondedEnergy = test_ljCos2.curEnergy['non_bonded', 0, 0] + \
            test_ljCos2.curEnergy['non_bonded', 0, 1] + \
            test_ljCos2.curEnergy['non_bonded', 1, 1]
        self.assertLessEqual(
            abs(totalEnergy - nonbondedEnergy), test_ljCos2.epsilon)

    def test_pressure(self):
        totalPressure = test_ljCos2.curPressure['total']
        nonbondedPressure = test_ljCos2.curPressure['non_bonded', 0, 0] + \
            test_ljCos2.curPressure['non_bonded', 0, 1] + \
            test_ljCos2.curPressure['non_bonded', 1, 1]
        self.assertLessEqual(
            abs(totalPressure - nonbondedPressure), test_ljCos2.epsilon)

    def test_relative_energy(self):
        relativeEnergyDiff = abs(
            (test_ljCos2.curEnergy['total'] - test_ljCos2.prevTotalEnergy) / test_ljCos2.prevTotalEnergy)
        self.assertLessEqual(relativeEnergyDiff, test_ljCos2.epsilon)

    def test_relative_pressure(self):
        relativePressureDiff = abs(
            (test_ljCos2.curPressure['total'] - test_ljCos2.prevTotalPressure) / test_ljCos2.prevTotalPressure)
        self.assertLessEqual(relativePressureDiff, test_ljCos2.epsilon)

    def test_max_force(self):
        # Checks the maximum deviation of Force is under epsilon
        maxdFx = 0.
        maxdFy = 0.
        maxdFz = 0.

        for i in range(0, test_ljCos2.numPart):
            dFx = abs(
                test_ljCos2.prevParticleData[i][5] - test_ljCos2.system.part[i].f[0])
            dFy = abs(
                test_ljCos2.prevParticleData[i][6] - test_ljCos2.system.part[i].f[1])
            dFz = abs(
                test_ljCos2.prevParticleData[i][7] - test_ljCos2.system.part[i].f[2])
            if(dFx > maxdFx):
                maxdFx = dFx
            if(dFy > maxdFy):
                maxdFy = dFy
            if(dFz > maxdFz):
                maxdFz = dFz

        self.assertLessEqual(maxdFx, 10 * test_ljCos2.epsilon)
        self.assertLessEqual(maxdFy, 10 * test_ljCos2.epsilon)
        self.assertLessEqual(maxdFz, 10 * test_ljCos2.epsilon)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    unittest.main()
