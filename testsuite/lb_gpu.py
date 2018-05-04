# Basic tests of the Lattice Boltzmann implementation       
#                                                           
# 1) check conservation of fluid mass                       
# 2) check conservation of total momentum                   
# 3) measure temperature of colloid and fluid  

from __future__ import print_function
import sys
import numpy as np
import unittest as ut
import espressomd
import espressomd.lb
from espressomd import *
from tests_common import abspath

@ut.skipIf(not espressomd.has_features(["LB_GPU","LENNARD_JONES"]) or espressomd.has_features("SHANCHEN"),
           "Features not available, skipping test!")
class TestLBGPU(ut.TestCase):
    box_l = 30.0
    system = espressomd.System(box_l=[box_l]*3)
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = range(n_nodes)

    def test(self):
        #setup parameters
        system = self.system
        int_steps = 5
        int_times = 1
        time_step = 0.005
        tau = 0.02
        agrid = 1.0
        dens = 0.85
        viscosity = 30.0
        friction = 2.0
        temp = 1.0
        gamma = 1.0
        skin = 0.4
        mom_prec = 1.e-2
        mass_prec_per_node = 4.e-8
        temp_confidence = 10

        if espressomd.has_features("ROTATION"):
            dof = 6.
        else:
            dof = 3.

        system.periodicity = [1, 1, 1]
        system.time_step = time_step
        system.cell_system.skin = skin

        #  Clear actors that might be left from prev tests
        if len(system.actors):
            del system.actors[0]
        system.part.clear()
        #import particle data
        data = np.genfromtxt(abspath("data/lb_system.data"))

        for particle in data:
            id = particle[0]
            typ = particle[1]
            pos = particle[3:6]
            f = particle[9:]
            v = particle[6:9]
            p=system.part.add(id=int(id), pos=pos, v=v, type=int(typ))
            if espressomd.has_features("ROTATION"):
                p.rotation=1,1,1

        n_col_part=len(system.part)

        system.time_step = time_step
        system.thermostat.set_langevin(kT=temp, gamma=gamma)
        system.integrator.run(100)
        #kill particle motion
        for i in range(n_col_part): 
            system.part[i].v=[0.0,0.0,0.0]
        system.thermostat.turn_off()
        
        lbf = lb.LBFluidGPU(visc=viscosity, dens=dens, agrid=agrid, tau=system.time_step, fric=friction)
        system.actors.add(lbf)
        system.thermostat.set_lb(kT=temp)
        #give particles a push
        for i in range(n_col_part): 
            system.part[i].v=system.part[i].v+[0.1,0.0,0.0]

        fluidmass = dens
        tot_mom = [0.0,0.0,0.0]
        for i in range(n_col_part):
            tot_mom=tot_mom + system.part[i].v

        system.integrator.run(100)

        max_dmass = 0.0
        max_dm = [0,0,0]
        avg_temp = 0.0
        avg_fluid_temp = 0.0
        
        #Integration
        for i in range(int_times):

            system.integrator.run(int_steps)

            fluidmass_sim=0.0
            fluid_temp_sim=0.0
            node_v_list = []
            node_dens_list = []
            for i in range(int(self.box_l/agrid)):
                for j in range(int(self.box_l/agrid)):
                    for k in range(int(self.box_l/agrid)):
                        node_v_list.append(lbf[i,j,k].velocity)
                        node_dens_list.append(lbf[i,j,k].density[0])

            #check mass conversation
            fluidmass_sim = sum(node_dens_list)/len(node_dens_list)

            dmass=abs(fluidmass_sim-fluidmass)#/len(node_dens_list)
            #self.assertTrue(dmass < mass_prec_per_node)
            if dmass > max_dmass:
                max_dmass=dmass     
            self.assertTrue(max_dmass < mass_prec_per_node, msg="fluid mass deviation too   high\ndeviation: {}   accepted deviation: {}".format(max_dmass, mass_prec_per_node))

            #check momentum conversation
            c_mom=system.analysis.analyze_linear_momentum()
            dm=abs(c_mom-tot_mom)
            #self.assertTrue(dm[0] <= mom_prec and dm[1] <= mom_prec and dm[2] <= mom_prec)
            for j in range(3):
                if dm[j] > max_dm[j]:
                    max_dm[j]=dm[j]
            self.assertTrue(max_dm[0] <= mom_prec and max_dm[1] <= mom_prec and max_dm[2] <= mom_prec, msg="momentum deviation too high\ndeviation: {}  accepted deviation: {}".format(max_dm, mom_prec))

            #check temp of particles
            e=system.analysis.energy()
            temp_particle = 2.0/dof*e["total"]/n_col_part
            avg_temp = avg_temp+temp_particle/int_times
            #check temp of fluid
            fluid_temp=0
            for j in range(len(node_dens_list)):
                fluid_temp = fluid_temp+(1.0/3.0)*(node_v_list[j][0]**2.0+node_v_list[j][1]**2.0+node_v_list[j][2]**2.0)*node_dens_list[j]*(self.box_l)**3/len(node_dens_list)**2
            avg_fluid_temp=avg_fluid_temp+fluid_temp/int_times



        temp_dev = (2.0/(n_col_part*3.0))**0.5
        temp_prec = temp_confidence*temp_dev/(int_times)**0.5
        
        print("maximal mass deviation: {}  accepted deviation: {}".format(max_dmass, mass_prec_per_node))
        print("maximal momentum deviation: {}  accepted deviation: {}".format(max_dm, mom_prec))
        print("average temperature: {}".format(avg_temp))
        print("average fluid temperature: {}".format(avg_fluid_temp))
        print("set temperature: {}".format(temp))
        print("maximally accepted deviation: {}".format(temp_prec))

        self.assertTrue(abs(avg_temp-temp) < temp_prec, msg="particle temperature deviation too high\ndeviation: {}  accepted deviation: {}".format(abs(avg_temp-temp), temp_prec))
        self.assertTrue(abs(avg_fluid_temp-temp) < temp_prec, msg="fluid temperature deviation too high\ndeviation: {}  accepted deviation: {}".format(abs(avg_fluid_temp-temp), temp_prec))


        
if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(TestLBGPU))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
