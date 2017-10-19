from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import FeneBond
from espressomd import polymer



@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeDistance(ut.TestCase):
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(1234)
    num_poly=4
    num_mono=20

    @classmethod
    def setUpClass(self):
        box_l = 100.0
        self.system.box_l = [box_l, box_l, box_l]
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        fene=FeneBond(k=30, d_r_max=2)
        self.system.bonded_inter.add(fene)
        polymer.create_polymer(N_P=self.num_poly,
                               bond_length=0.9,
                               MPC=self.num_mono,
                               bond=fene)

    # python version of the espresso core function,
    # this function assumes the polymer is entirely contained in a sim box
    def calc_re(self):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        dist = np.fabs(self.system.part[head_id].pos - self.system.part[tail_id].pos)
        # check smaller distances via PBC
        dist = np.where(dist > 0.5 * self.system.box_l, self.system.box_l - dist, dist)
        dist = np.sum(dist**2, axis=-1)
        return np.mean(np.sqrt(dist)),np.std(np.sqrt(dist)), np.mean(dist), np.std(dist)

    # python version of the espresso core function,
    # this function assumes the polymer is entirely contained in a sim box
    def calc_rg(self):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        rg2=[]
        for p in range(self.num_poly):
            rg2.append(np.var(self.system.part[head_id[p]:tail_id[p]+1].pos, axis=0))
        rg2 = np.array(rg2)
        rg2 = np.sum(rg2, axis=1)
        return np.mean(np.sqrt(rg2)), np.std(np.sqrt(rg2)), np.mean(rg2), np.std(rg2)

    # python version of the espresso core function,
    # this function assumes the polymer is entirely contained in a sim box
    def calc_rh(self):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        rh=[]
        for p in range(self.num_poly):
            r = np.array(self.system.part[head_id[p]:tail_id[p]+1].pos)
            # this generates indices for all i<j combinations
            ij = np.triu_indices(len(r), k=1)
            r_ij = np.fabs(r[ij[0]] - r[ij[1]])
            # check smaller distances via PBC
            r_ij = np.where(r_ij > 0.5 * self.system.box_l, self.system.box_l-r_ij, r_ij)
            dist = np.sqrt(np.sum(r_ij**2, axis=1))
            rh.append(self.num_mono*self.num_mono*0.5/(np.sum(1./dist)))
            # the other way do it, with the proper prefactor of N(N-1)
            #rh.append(np.mean(1./dist))
        rh = np.array(rh)
        return np.mean(rh), np.std(rh)


    def test_radii(self):
        # compare calc_re()
        core_re=self.system.analysis.calc_re(chain_start=0,
                                              number_of_chains=self.num_poly,
                                              chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_re, self.calc_re()))

        # compare calc_rg()
        core_rg=self.system.analysis.calc_rg(chain_start=0,
                                              number_of_chains=self.num_poly,
                                              chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_rg, self.calc_rg()))

        # compare calc_rh()
        core_rh=self.system.analysis.calc_rh(chain_start=0,
                                              number_of_chains=self.num_poly,
                                              chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_rh, self.calc_rh()))




if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
