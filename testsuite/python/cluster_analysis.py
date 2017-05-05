#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond
from espressomd.pair_criteria import *
from espressomd.cluster_analysis import ClusterStructure


class ClusterAnalysis(ut.TestCase):
    """Tests the cluster analysis"""

    es = espressomd.System()

    f = FeneBond(k=1, d_r_max=5)
    es.bonded_inter.add(f)
    
    # Firt cluster
    es.part.add(id=0,pos=(0,0,0))
    es.part.add(id=1,pos=(0.91,0,0),bonds=((0,0),))
    es.part.add(id=2,pos=(0,0.2,0))
    es.part.add(id=3,pos=(0,0.1,0))

    # 2nd cluster
    es.part.add(id=4,pos=(0.5,0.5,0.5))
    es.part.add(id=5,pos=(0.55,0.5,0.5))
    cs=ClusterStructure()

    def test_00_fails_without_criterion_set(self):
        self.assertRaises(self.cs.run_for_all_pairs())
    
    def test_set_criterion(self):
        # Test setters/getters for criteria
        dc=DistanceCriterion(cut_off=0.11)
        self.cs.set_params(pair_criterion=dc)
        # Do we get back the right criterion
        dc_ret =self.cs.get_params()["pair_criterion"]
        # Note: This work around the fact that the script interface does not 
        # yet assign the correct derived class when returning an object
        self.assertTrue(dc_ret.name() == "PairCriteria::DistanceCriterion")
        self.assertTrue(abs(dc_ret.get_params()["cut_off"]-0.11)<=1E-8)

        # Is the cluster structure empty before being used

    def test_analysis_for_all_pairs(self):
        # Run cluster analysis 
        self.cs.set_params(pair_criterion=DistanceCriterion(cut_off=0.12))
        self.cs.run_for_all_pairs()

        # Number of clusters
        self.assertTrue(len(self.cs.clusters)==2)
        cids=self.cs.cluster_ids()
        print(cids)
        print(self.cs.clusters)
        print(self.cs.clusters[cids[0]])

        # Sizes of individual clusters
        l2=self.cs.clusters[cids[1]].size()
        l1=len(self.cs.clusters[cids[0]].particle_ids())
        
        # Clusters should contain 2 and 4 particles
        self.assertTrue(min(l1,l2)==2)
        self.assertTrue(max(l1,l2)==4)

        # Verify particle ids
        smaller_cluster=None
        bigger_cluster=None
        if l1<l2:
            smaller_cluster=self.cs.clusters[cids[0]]
            bigger_cluster=self.cs.clusters[cids[1]]
        else:
            smaller_cluster=self.cs.clusters[cids[1]]
            bigger_cluster=self.cs.clusters[cids[0]]

        self.assertTrue(bigger_cluster.particle_ids()==[0,1,2,3])
        self.assertTrue(smaller_cluster.particle_ids()==[4,5])

        # Test obtaining a ParticleSlice for a cluster
        pids=bigger_cluster.particle_ids()
        particles=bigger_cluster.particles()
        # Do the number of entries match
        self.assertTrue(len(pids)==len(particles.id_selection))
        
        # Compare ids of particles in the slice
        self.assertTrue(all(particles.id_selection)==all(pids))


        # Test iteration over clusters
        visited_sizes=[]
        for c in self.cs.clusters:
            visited_sizes.append(c[1].size())
        visited_sizes=sorted(visited_sizes)
        self.assertTrue(visited_sizes==[2,4])
        


        
        
          
    
    def test_analysis_for_bonded_particles(self):
        # Run cluster analysis 
        self.cs.set_params(pair_criterion=BondCriterion(bond_type=0))
        self.cs.run_for_bonded_particles()

        # There should be one cluster containing particles 0 and 1
        self.assertTrue(len(self.cs.clusters)==1)
        self.assertTrue(self.cs.clusters[self.cs.cluster_ids()[0]].particle_ids()==[0,1])
        
        # Check particle to cluster id mapping, once by ParticleHandle, once by id
        self.assertTrue(self.cs.cid_for_particle(self.es.part[0])==self.cs.cluster_ids()[0])
        self.assertTrue(self.cs.cid_for_particle(1)==self.cs.cluster_ids()[0])

        

if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
