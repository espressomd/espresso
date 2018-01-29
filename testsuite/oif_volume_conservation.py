import espressomd


import numpy as np
import os
import sys

import unittest as ut
from tests_common import abspath


@ut.skipIf(not espressomd.has_features("MEMBRANE_COLLISION", "OIF_LOCAL_FORCES", "OIF_GLOBAL_FORCES"), "OIF featues not compiled in.")
class OifVolumeConservation(ut.TestCase):

    """Loads a soft elastic sphere via object_in_fluid, stretches it and checks resotration of original volume due to elastic forces."""

    def test(self):
        import object_in_fluid as oif

        system = espressomd.System()
        self.assertEqual(system.max_oif_objects,0)
        system.time_step = 0.1
        system.box_l = [10.0, 10.0, 10.0]
        system.cell_system.skin = 0.5
        system.thermostat.set_langevin(kT=0, gamma=0.7)

        # creating the template for OIF object
        cell_type = oif.OifCellType(nodes_file=abspath("data/sphere393nodes.dat"), triangles_file=abspath(
            "data/sphere393triangles.dat"), system=system, ks=1.0, kb=1.0, kal=1.0, kag=0.1, kv=0.1, check_orientation=False, resize=(3.0, 3.0, 3.0))

        # creating the OIF object
        cell0 = oif.OifCell(
            cell_type=cell_type, part_type=0, origin=[5.0, 5.0, 5.0])
        self.assertEqual(system.max_oif_objects,1)
        # cell0.output_vtk_pos_folded(file_name="cell0_0.vtk")

        # fluid

        diameter_init = cell0.diameter()
        print("initial diameter = " + str(diameter_init))

        # OIF object is being stretched by factor 1.5
        maxCycle = 500
        system.part[:].pos = (system.part[:].pos - 5) * 1.5 + 5

        diameter_stretched = cell0.diameter()
        print("stretched diameter = " + str(diameter_stretched))

        # main integration loop
        # OIF object is let to relax into relaxed shape of the sphere
        for i in range(3):
            system.integrator.run(steps=300)
            diameter_final = cell0.diameter()
            print("final diameter = " + str(diameter_final))
            self.assertAlmostEqual(
                diameter_final / diameter_init - 1, 0, delta=0.005)
        


if __name__ == "__main__":
    # print("Features: ", espressomd.features())
    ut.main()
