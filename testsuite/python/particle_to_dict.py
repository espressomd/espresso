import unittest as ut
import unittest_decorators as utx
import numpy as np
import espressomd

class ParticleDictionaryTest(ut.TestCase):
    box_l = 30.
    system = espressomd.System(box_l = 3 * [box_l])

    def test(self):
        p = self.system.part.add(pos = np.random.uniform(size = (10,3)) * self.box_l)
        pp = str(p)
        pdict = p.to_dict()
        p.remove()
        self.system.part.add(pdict)
        self.assertEqual(str(self.system.part.select()), pp)

if __name__ == "__main__":
    ut.main()
