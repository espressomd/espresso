import numpy as np
import unittest as ut

import espressomd
import espressomd.electrostatics 
import espressomd.interactions
from espressomd import drude_helpers

class Drude(ut.TestCase):

    @ut.skipIf(not espressomd.has_features("P3M", "THOLE", "LANGEVIN_PER_PARTICLE"), "Test needs P3M, THOLE and LANGEVIN_PER_PARTICLE")
    def test(self):
        """
        Sets up a BMIM PF6 pair separated in y-direction with fixed cores.
        Adds the Drude particles and related features (intramolecular exclusion bonds, Thole screening)
        via helper functions.
        Calculates the induced dipole moment and the diagonals of the polarization tensor
        and compares against reference results, which where reproduced with LAMMPS.

        """
        
        box_l = 50
        system=espressomd.System(box_l = [box_l, box_l, box_l])

        system.seed = system.cell_system.get_state()['n_nodes'] * [12]
        np.random.seed(12)

        #Reference Results, reproduced with LAMMPS
        #Dipole Moments
        ref_mu0_pf6  = [ 0.00177594, 0.16480996, -0.01605161] 
        ref_mu0_c1   = [ 0.00076652, 0.15238767,  0.00135291] 
        ref_mu0_c2   = [-0.00020222, 0.11084197,  0.00135842] 
        ref_mu0_c3   = [ 0.00059177, 0.23949626, -0.05238468] 
        ref_mu0_bmim = [ 0.00115606, 0.5027259 , -0.04967335] 

        #Polarisation Tensor diagonals
        ref_pol_pf6  = [4.5535698335873445, 4.7558611769477697, 4.5546580162000554]
        ref_pol_bmim = [13.126868394164262, 14.392582501485913, 16.824150151623762]

        #TIMESTEP
        fs_to_md_time = 1.0e-2
        time_step_fs = 0.5
        time_step_ns = time_step_fs*1e-6
        dt = time_step_fs * fs_to_md_time

        #COM TEMPERATURE
        #Global thermostat temperature, for com and langevin. 
        #LangevinPerParticle temperature is set to 0 for drude and core to properly account for com forces. 
        #Like that, langevin thermostat can still be used for non-drude particles
        SI_temperature = 300.0
        gamma_com = 1.0
        kb_kjmol = 0.0083145
        temperature_com = SI_temperature*kb_kjmol

        #COULOMB PREFACTOR (elementary charge)^2 / (4*pi*epsilon_0) in Angstrom * kJ/mol
        coulomb_prefactor = 1.67101e5*kb_kjmol

        #POLARIZATION
        #polarization  = 1.0 #In (Angstrom^3)_CGS 
        #alpha_SI = 4*Pi*eps_0 alpha_CGS; 4*Pi*epsilon_0*Angstrom^3/((elementary charge)^2*Angstrom^2*N_A/kJ)
        conv_pol_CGS_SI = 7.197586e-4
        #alpha = conv_pol_CGS_SI*args.polarization

        #DRUDE/TOTAL MASS
        #lamoureux03 used values 0.1-0.8 g/mol for drude mass
        mass_drude = 0.8
        mass_tot = 100.0
        mass_core = mass_tot-mass_drude
        mass_red_drude = mass_drude*mass_core/mass_tot

        #SPRING CONSTANT DRUDE
        #Used 1000kcal/mol/A^2 from lamoureux03a table 1 p 3031
        k_drude = 4184.0;#in kJ/mol/A^2
        T_spring = 2.0*np.pi*np.sqrt(mass_drude/k_drude)
        #T_spring_fs = T_spring/fs_to_md_time
        #Period of free oscillation: T_spring = 2Pi/w; w = sqrt(k_d/m_d)

        #TEMP DRUDE
        #Used T* = 1K from lamoureux03a p 3031 (2) 'Cold drude oscillators regime'
        SI_temperature_drude = 1.0
        temperature_drude = SI_temperature_drude*kb_kjmol

        #GAMMA DRUDE
        #Thermostat relaxation time should be similar to T_spring
        gamma_drude = mass_red_drude/T_spring

        system.cell_system.skin=0.4
        system.time_step=dt

        #Forcefield
        types           = {"PF6":          0, "BMIM_C1":           1, "BMIM_C2":         2, "BMIM_C3":          3, "BMIM_COM":  4, "PF6_D": 5, "BMIM_C1_D": 6, "BMIM_C2_D": 7, "BMIM_C3_D": 8}
        charges         = {"PF6":      -0.78, "BMIM_C1":      0.4374, "BMIM_C2":    0.1578, "BMIM_C3":     0.1848, "BMIM_COM":  0}
        polarizations   = {"PF6":      4.653, "BMIM_C1":       5.693, "BMIM_C2":     2.103, "BMIM_C3":      7.409} 
        masses          = {"PF6":     144.96, "BMIM_C1":       67.07, "BMIM_C2":     15.04, "BMIM_C3":      57.12, "BMIM_COM":  0}
        masses["BMIM_COM"] = masses["BMIM_C1"] +  masses["BMIM_C2"] +  masses["BMIM_C3"]

        box_center = 0.5*np.array(3*[box_l])
        system.min_global_cut = 3.5

        #Place Particles
        dmol = 5.0

        #Test Anion
        pos_pf6 = box_center + np.array([0,dmol,0])
        system.part.add(id=0, type=types["PF6"],  pos = pos_pf6, q=charges["PF6"], mass = masses["PF6"], fix = [1,1,1])

        pos_com = box_center - np.array([0,dmol,0])
        system.part.add(id = 2, type = types["BMIM_C1"], pos = pos_com + [0, -0.527, 1.365], q = charges["BMIM_C1"], mass=masses["BMIM_C1"], fix = [1,1,1])
        system.part.add(id = 4, type = types["BMIM_C2"], pos = pos_com + [0, 1.641, 2.987], q = charges["BMIM_C2"], mass=masses["BMIM_C2"], fix = [1,1,1])
        system.part.add(id = 6, type = types["BMIM_C3"], pos = pos_com + [0, 0.187, -2.389], q = charges["BMIM_C3"], mass=masses["BMIM_C3"], fix = [1,1,1])

        system.thermostat.set_langevin(kT=temperature_com, gamma=gamma_com)

        p3m=espressomd.electrostatics.P3M(prefactor=coulomb_prefactor, accuracy=1e-4, mesh = [18,18,18], cao = 5)

        system.actors.add(p3m)

        #Drude related Bonds

        thermalized_dist_bond = espressomd.interactions.ThermalizedBond(temp_com = temperature_com, gamma_com = gamma_com, temp_distance = temperature_drude, gamma_distance = gamma_drude, r_cut = 1.0)
        harmonic_bond = espressomd.interactions.HarmonicBond(k = k_drude, r_0 = 0.0, r_cut = 1.0)
        system.bonded_inter.add(thermalized_dist_bond)
        system.bonded_inter.add(harmonic_bond)

        drude_helpers.add_drude_particle_to_core(system, harmonic_bond, thermalized_dist_bond, system.part[0], 1, types["PF6_D"], polarizations["PF6"], mass_drude, coulomb_prefactor, 2.0)
        drude_helpers.add_drude_particle_to_core(system, harmonic_bond, thermalized_dist_bond, system.part[2], 3, types["BMIM_C1_D"], polarizations["BMIM_C1"], mass_drude, coulomb_prefactor, 2.0)
        drude_helpers.add_drude_particle_to_core(system, harmonic_bond, thermalized_dist_bond, system.part[4], 5, types["BMIM_C2_D"], polarizations["BMIM_C2"], mass_drude, coulomb_prefactor, 2.0)
        drude_helpers.add_drude_particle_to_core(system, harmonic_bond, thermalized_dist_bond, system.part[6], 7, types["BMIM_C3_D"], polarizations["BMIM_C3"], mass_drude, coulomb_prefactor, 2.0)

        #Setup and add Drude-Core SR exclusion bonds
        drude_helpers.setup_and_add_drude_exclusion_bonds(system)

        #Setup intramol SR exclusion bonds once
        drude_helpers.setup_intramol_exclusion_bonds(system, [6,7,8],[1,2,3], [charges["BMIM_C1"], charges["BMIM_C2"], charges["BMIM_C3"]])

        #Add bonds per molecule
        drude_helpers.add_intramol_exclusion_bonds(system, [3,5,7], [2,4,6])

        #Thole
        drude_helpers.add_all_thole(system)

        def dipole_moment(id_core, id_drude):
            pc = system.part[id_core]
            pd = system.part[id_drude]
            v = pd.pos - pc.pos
            return pd.q * v

        def measure_dipole_moments():
            dm_pf6 = []
            dm_C1 = []
            dm_C2 = []
            dm_C3 = []
            
            system.integrator.run(115)
            
            for i in range(100): 
                system.integrator.run(1)
                
                dm_pf6.append(dipole_moment(0,1))
                dm_C1.append(dipole_moment(2,3))
                dm_C2.append(dipole_moment(4,5))
                dm_C3.append(dipole_moment(6,7))
            
            dm_pf6_m = np.mean(dm_pf6, axis = 0)
            dm_C1_m = np.mean(dm_C1, axis = 0)
            dm_C2_m = np.mean(dm_C2, axis = 0)
            dm_C3_m = np.mean(dm_C3, axis = 0)
            dm_sum_bmim = dm_C1_m + dm_C2_m + dm_C3_m

            res = dm_pf6_m, dm_C1_m, dm_C2_m, dm_C3_m, dm_sum_bmim
            return res


        def setElectricField(E):
            E=np.array(E)
            for p in system.part:
                p.ext_force = p.q*E

        def calc_pol(mu0, muE, E):
            pol = (muE-mu0)/E/conv_pol_CGS_SI
            return pol

        def measure_pol(Es, dim):
            E = [0.0,0.0,0.0]
            E[dim] = Es

            setElectricField(E)
            mux_pf6, mux_c1, mux_c2, mux_c3, mux_bmim = measure_dipole_moments()

            return calc_pol(mu0_pf6[dim], mux_pf6[dim], Es), calc_pol(mu0_bmim[dim], mux_bmim[dim], Es)

        mu0_pf6, mu0_c1, mu0_c2, mu0_c3, mu0_bmim = measure_dipole_moments() 
        eA_to_Debye=4.8032047
        atol = 1e-2
        rtol = 1e-2

        np.testing.assert_allclose(ref_mu0_pf6,  eA_to_Debye * mu0_pf6, atol = atol, rtol = rtol) 
        np.testing.assert_allclose(ref_mu0_c1,   eA_to_Debye * mu0_c1, atol = atol, rtol = rtol) 
        np.testing.assert_allclose(ref_mu0_c2,   eA_to_Debye * mu0_c2, atol = atol, rtol = rtol) 
        np.testing.assert_allclose(ref_mu0_c3,   eA_to_Debye * mu0_c3, atol = atol, rtol = rtol) 
        np.testing.assert_allclose(ref_mu0_bmim, eA_to_Debye * mu0_bmim, atol = atol, rtol = rtol) 

        pol_pf6 = []
        pol_bmim = []
        Efield=96.48536  # = 1 V/A in kJ / (Avogadro Number) / Angstrom / elementary charge

        res = measure_pol(Efield, 0)
        pol_pf6.append(res[0])
        pol_bmim.append(res[1])

        res = measure_pol(Efield, 1)
        pol_pf6.append(res[0])
        pol_bmim.append(res[1])

        res = measure_pol(Efield, 2)
        pol_pf6.append(res[0])
        pol_bmim.append(res[1])

        np.testing.assert_allclose(ref_pol_pf6, pol_pf6, atol = atol, rtol = rtol) 
        np.testing.assert_allclose(ref_pol_bmim, pol_bmim, atol = atol, rtol = rtol ) 

if __name__ == "__main__":
    ut.main()
