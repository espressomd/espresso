from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import FeneBond
from espressomd import polymer

@ut.skipIf(not espressomd.has_features("LENNARD_JONES"), "Skipped because LENNARD_JONES turned off.")
class AnalyzeChain(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1234)
    num_poly=2
    num_mono=5

    @classmethod
    def setUpClass(self):
        box_l = 20.0
        # start with a small bo
        self.system.box_l = np.array([box_l, box_l, box_l])
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        fene=FeneBond(k=30, d_r_max=2)
        self.system.bonded_inter.add(fene)
        polymer.create_polymer(N_P=self.num_poly,
                               bond_length=0.9,
                               MPC=self.num_mono,
                               bond=fene)
        # bring two polymers to opposite corners:
        # far in centre cell, but mirror images are close
        head_id = 0
        tail_id = head_id+self.num_mono
        cm=np.mean(self.system.part[head_id:tail_id].pos, axis=0)
        self.system.part[head_id:tail_id].pos = self.system.part[head_id:tail_id].pos - cm + self.system.box_l
        head_id = self.num_mono+1
        tail_id = head_id+self.num_mono
        cm=np.mean(self.system.part[head_id:tail_id].pos, axis=0)
        self.system.part[head_id:tail_id].pos -= cm

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_re(self):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        dist = self.system.part[head_id].pos - self.system.part[tail_id].pos
        dist = np.sum(dist**2, axis=-1)
        return np.mean(np.sqrt(dist)),np.std(np.sqrt(dist)), np.mean(dist), np.std(dist)

    # python version of the espresso core function,
    # does not check mirror distances
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
    # does not check mirror distances
    def calc_rh(self):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        rh=[]
        for p in range(self.num_poly):
            r = np.array(self.system.part[head_id[p]:tail_id[p]+1].pos)
            # this generates indices for all i<j combinations
            ij = np.triu_indices(len(r), k=1)
            r_ij = r[ij[0]] - r[ij[1]]
            dist = np.sqrt(np.sum(r_ij**2, axis=1))
            #rh.append(self.num_mono*self.num_mono*0.5/(np.sum(1./dist)))
            # the other way do it, with the proper prefactor of N(N-1)
            rh.append(1./np.mean(1./dist))
        rh = np.array(rh)
        return np.mean(rh), np.std(rh)

    # python version of the espresso core function,
    # does not check mirror distances
    def calc_chain_rdf(self,bins):
        head_id=np.arange(0,self.num_poly*self.num_mono,self.num_mono)
        tail_id=head_id+self.num_mono-1
        data_cm=[]
        data_min=[]
        data_mono=[]
        i=np.arange(self.num_mono*self.num_mono)//self.num_mono
        j=np.arange(self.num_mono*self.num_mono)%self.num_mono
        for p1 in range(self.num_poly):
            for p2 in range(p1):
                r1 = np.array(self.system.part[head_id[p1]:tail_id[p1]+1].pos)
                r2 = np.array(self.system.part[head_id[p2]:tail_id[p2]+1].pos)
                r_ij = r1[i] - r2[j]
                dist = np.sqrt(np.sum(r_ij**2, axis=1))
                data_min.append(np.min(dist))
                dist_cm = np.sqrt(np.sum((np.mean(r1,axis=0)-np.mean(r2,axis=0))**2))
                data_cm.append(dist_cm)
                data_mono.append(dist)
        hist_mono = np.histogram(np.ravel(data_mono), bins=bins, density=False)[0]
        hist_cm =np.histogram(data_cm, bins=bins, density=False)[0]
        hist_min =np.histogram(data_min, bins=bins, density=False)[0]
        return hist_mono, hist_cm, hist_min

    # test core results versus python variants (no PBC)
    def test_radii(self):
        # increase PBC for remove mirror images
        old_pos = self.system.part[:].pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.system.part[:].pos = old_pos
        # compare calc_re()
        core_re = self.system.analysis.calc_re(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_re, self.calc_re()))
        # compare calc_rg()
        core_rg = self.system.analysis.calc_rg(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_rg, self.calc_rg()))
        # compare calc_rh()
        core_rh = self.system.analysis.calc_rh(chain_start=0,
                                               number_of_chains=self.num_poly,
                                               chain_length=self.num_mono)
        self.assertTrue( np.allclose(core_rh, self.calc_rh()))
        # restore PBC
        self.system.box_l = self.system.box_l / 2.
        self.system.part[:].pos = old_pos

    # test core results versus python variants (no PBC)
    def test_chain_rdf(self):
        # increase PBC for remove mirror images
        old_pos = self.system.part[:].pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.system.part[:].pos = old_pos
        r_min = 0.0
        r_max = 100.0
        r_bins = 10
        bin_width = (r_max-r_min)/r_bins
        bins = np.arange(r_min, r_max+bin_width, bin_width)
        core_rdf = self.system.analysis.rdf_chain(r_min=r_min,
                                                 r_max=r_max,
                                                 r_bins=r_bins,
                                                 chain_start=0,
                                                 number_of_chains=self.num_poly,
                                                 chain_length=self.num_mono)
        rdf=(self.calc_chain_rdf(bins))
        bin_volume = 4./3.*np.pi*(bins[1:]**3 -  bins[:-1]**3)
        box_volume = self.system.box_l[0]*self.system.box_l[1]*self.system.box_l[2]
        # number pairs between in all particles
        num_pair_part = 0.5*(self.num_mono*self.num_poly)*(self.num_mono*self.num_poly-1)
        # number of polymer pairs
        num_pair_poly = 0.5*(self.num_poly)*(self.num_poly-1)
        # number of pairs between monomers or different polymers
        num_pair_mono = 0.5*(self.num_mono*self.num_mono)*(self.num_poly-1)*(self.num_poly)
        # bins
        self.assertTrue( np.allclose(core_rdf[:,0], (bins[1:]+bins[:-1])*0.5))
        # monomer rdf
        self.assertTrue( np.allclose(core_rdf[:,1]*bin_volume*num_pair_mono/box_volume, rdf[0]))
        # cm rdf
        self.assertTrue( np.allclose(core_rdf[:,2]*bin_volume*num_pair_poly/box_volume, rdf[1]))
        # min rdf
        self.assertTrue( np.allclose(core_rdf[:,3]*bin_volume*num_pair_poly/box_volume, rdf[2]))
        # restore PBC
        self.system.box_l = self.system.box_l / 2.
        self.system.part[:].pos = old_pos

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
