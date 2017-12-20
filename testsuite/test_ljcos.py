# This is a unit-test for Lennard Jones Cosine non-bonded interaction

from __future__ import print_function
import unittest
import espressomd
import numpy as np


@unittest.skipIf(not espressomd.has_features(["LJCOS"]),
                 "Features not available, skipping test!")
class test_ljCos(unittest.TestCase):

    test_file = 'data/test_ljcos.dat'

    # Test Parameters
    epsilon = 1E-4
    numPart = 512

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
        file = open(cls.test_file, 'r')
        for line in file:
            cls.data.append(line.split())
        cls.prevTotalEnergy = float(cls.data[3][0])
        cls.prevTotalPressure = float(cls.data[3][1])
        cls.system.box_l = map(
            float, [cls.data[5][1], cls.data[5][2], cls.data[5][3]])
        for i in range(9, len(cls.data)):
            cls.system.part.add(id=int(cls.data[i][0]), type=int(cls.data[i][1]),
                                pos=map(float, [cls.data[i][2], cls.data[i][3], cls.data[i][4]]))
            cls.prevParticleData.append(map(float, cls.data[i]))

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
        cls.system.non_bonded_inter[0, 0].lennard_jones_cos.set_params(
            epsilon=1., sigma=1.,
            cutoff=1.12246, offset=0.0)
        cls.system.non_bonded_inter[0, 1].lennard_jones_cos.set_params(
            epsilon=1.3, sigma=0.5,
            cutoff=2.0, offset=0.0)
        cls.system.non_bonded_inter[1, 1].lennard_jones_cos.set_params(
            epsilon=2.2, sigma=1.0,
            cutoff=1.12246, offset=0.5)
        cls.system.integrator.run(0)

    @classmethod
    def setUpClass(cls):
        cls.read()
        cls.time_translation_dt()
        # After time translation
        cls.curEnergy = cls.system.analysis.energy()
        cls.curPressure = cls.system.analysis.pressure()

    def test_energy(self):
        totalEnergy = test_ljCos.curEnergy['total']
        nonbondedEnergy = test_ljCos.curEnergy['non_bonded', 0, 0] + \
            test_ljCos.curEnergy['non_bonded', 0, 1] + \
            test_ljCos.curEnergy['non_bonded', 1, 1]
        self.assertLessEqual(
            abs(totalEnergy - nonbondedEnergy), test_ljCos.epsilon)

    def test_pressure(self):
        totalPressure = test_ljCos.curPressure['total']
        nonbondedPressure = test_ljCos.curPressure['non_bonded', 0, 0] + \
            test_ljCos.curPressure['non_bonded', 0, 1] + \
            test_ljCos.curPressure['non_bonded', 1, 1]
        self.assertLessEqual(
            abs(totalPressure - nonbondedPressure), test_ljCos.epsilon)

    def test_relative_energy(self):
        relativeEnergyDiff = abs(
            (test_ljCos.curEnergy['total'] - test_ljCos.prevTotalEnergy) / test_ljCos.prevTotalEnergy)
        self.assertLessEqual(relativeEnergyDiff, test_ljCos.epsilon)

    def test_relative_pressure(self):
        relativePressureDiff = abs(
            (test_ljCos.curPressure['total'] - test_ljCos.prevTotalPressure) / test_ljCos.prevTotalPressure)
        self.assertLessEqual(relativePressureDiff, test_ljCos.epsilon)

    def test_max_force(self):
        # Checks the maximum deviation of Force is under epsilon
        maxdFx = 0.
        maxdFy = 0.
        maxdFz = 0.

        for i in range(0, test_ljCos.numPart):
            dFx = abs(
                test_ljCos.prevParticleData[i][5] - test_ljCos.system.part[i].f[0])
            dFy = abs(
                test_ljCos.prevParticleData[i][6] - test_ljCos.system.part[i].f[1])
            dFz = abs(
                test_ljCos.prevParticleData[i][7] - test_ljCos.system.part[i].f[2])
            if(dFx > maxdFx):
                maxdFx = dFx
            if(dFy > maxdFy):
                maxdFy = dFy
            if(dFz > maxdFz):
                maxdFz = dFz

        self.assertLessEqual(maxdFx, 10 * test_ljCos.epsilon)
        self.assertLessEqual(maxdFy, 10 * test_ljCos.epsilon)
        self.assertLessEqual(maxdFz, 10 * test_ljCos.epsilon)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    unittest.main()
