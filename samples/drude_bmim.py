from __future__ import print_function
import sys
import time
import espressomd
from espressomd.electrostatics import P3M, P3M_GPU
from espressomd.interactions import DrudeBond
import os
import numpy as np
import argparse
from espressomd.visualization_opengl import *
from threading import Thread
from time import sleep
from drude_functions import *

parser = argparse.ArgumentParser(description='Drude LJ liquid')
parser.add_argument("--epsilon_r", nargs='?', default=1.0,type=float)
parser.add_argument("--mass_drude", nargs='?', default=0.8,type=float)
parser.add_argument("--walltime", nargs='?', default=1.0,type=float)
parser.add_argument("--drude", dest ="drude",action="store_true")
parser.add_argument("--no-drude", dest ="drude",action="store_false")
parser.add_argument("--thole", dest ="thole",action="store_true")
parser.add_argument("--no-thole", dest ="thole",action="store_false")
parser.add_argument("--intra_ex", dest ="intra_ex",action="store_true")
parser.add_argument("--no-intra_ex", dest ="intra_ex",action="store_false")
parser.add_argument("--visual", dest ="visu",action="store_true")
parser.add_argument("--gpu", dest ="gpup3m",action="store_true")
parser.add_argument("--timing", dest ="timing",action="store_true")
parser.set_defaults(drude=True, thole=True, intra_ex=True, visu=False, gpup3m=False, timing=False)
args = parser.parse_args()

print(args)

S=espressomd.System()

if args.visu:
    a = 0.2
    d_scale= 0.53
    c_ani = [1,0,0,1]
    c_dru = [0,0,1,a]
    c_com = [0,0,0,1]
    c_cat = [0,1,0,1]
    visualizer = openGLLive(S, 
    background_color=[1,1,1], 
    drag_enabled= True, 
    ext_force_arrows =True,
    drag_force = 10, 
    draw_bonds=False,
    quality_particles=32,
    particle_coloring = 'type',
    particle_type_colors = [c_ani,c_cat,c_cat,c_cat,c_com,c_dru,c_dru,c_dru,c_dru],
    particle_sizes = [0.5*5.06, 0.5*4.38, 0.5*3.41, 0.5*5.04, 0.1, d_scale*5.06, d_scale*4.38, d_scale*3.41, d_scale*5.04])

tpath=""
if args.thole:
    tpath = "_thole"

if args.drude:
    path = "/work/konrad/drude_python/BMIMPF6_bulk/"
else:
    path = "/work/konrad/drude_python/BMIMPF6_bulk_nodrude/"

if not os.path.exists(path):
    os.mkdir(path)

print("Writing to",path)

#Drude related parameter
#----------------------

#TIMESTEP
fs_to_md_time = 1.0e-2
time_step_fs = 2.0
time_step_ns = time_step_fs*1e-6
dt = time_step_fs * fs_to_md_time

#COM TEMPERATURE
#Global thermostat temperature, for com and langevin. 
#LangevinPerParticle temperature is set to 0 for drude and core to properly account for com forces. 
#Like that, langevin thermostat can still be used for non-drude particles
SI_temperature = 353.0
gamma_com = 1.0
kb_kjmol = 0.0083145
temperature_com = SI_temperature*kb_kjmol

#BJERRUM LENGTH
l_b  = 1.67101e-5/SI_temperature*1e10/args.epsilon_r
print("l_bjerrum:",l_b)

#POLARIZATION
#polarization  = 1.0 #In (Angstrom^3)_CGS 
#alpha_SI = 4*Pi*eps_0 alpha_CGS; 4*Pi*epsilon_0*Angstrom^3/((elementary charge)^2*Angstrom^2*N_A/kJ)
conv_pol_CGS_SI = 7.197586e-4
#alpha = conv_pol_CGS_SI*args.polarization

#DRUDE/TOTAL MASS
#lamoureux03 used values 0.1-0.8 g/mol for drude mass
mass_tot = 100.0
mass_core = mass_tot-args.mass_drude
mass_red_drude = args.mass_drude*mass_core/mass_tot

#SPRING CONSTANT DRUDE
#Used 1000kcal/mol/A^2 from lamoureux03a table 1 p 3031
k_drude = 4184.0;#in kJ/mol/A^2
T_spring = 2.0*np.pi*np.sqrt(args.mass_drude/k_drude)
#T_spring_fs = T_spring/fs_to_md_time
#Period of free oscillation: T_spring = 2Pi/w; w = sqrt(k_d/m_d)

#TEMP DRUDE
#Used T* = 1K from lamoureux03a p 3031 (2) 'Cold drude oscillators regime'
SI_temperature_drude = 1.0
temperature_drude = SI_temperature_drude*kb_kjmol

#GAMMA DRUDE
#Thermostat relaxation time should be similar to T_spring
gamma_drude = mass_red_drude/T_spring
print("Gamma_drude:",gamma_drude)

#System
deg_free = 3

S.cell_system.skin=0.4
S.time_step=dt


#Forcefield
types           = {"PF6":          0, "BMIM_C1":           1, "BMIM_C2":         2, "BMIM_C3":          3, "BMIM_COM":  4, "PF6_D": 5, "BMIM_C1_D": 6, "BMIM_C2_D": 7, "BMIM_C3_D": 8}
inv_types = {v: k for k, v in types.iteritems()}
charges         = {"PF6":      -0.78, "BMIM_C1":      0.4374, "BMIM_C2":    0.1578, "BMIM_C3":     0.1848, "BMIM_COM":  0}
#charges         = {"PF6":      -1.0, "BMIM_C1":      0.5608, "BMIM_C2":    0.2023, "BMIM_C3":     0.2369, "BMIM_COM":  0}
#polarizations   = {"PF6":     5.1825, "BMIM_C1":      5.8585, "BMIM_C2":    2.0439, "BMIM_C3":      7.3422} #Hannes
polarizations   = {"PF6":      4.653, "BMIM_C1":       5.693, "BMIM_C2":     2.103, "BMIM_C3":      7.409} #Frank
masses          = {"PF6":     144.96, "BMIM_C1":       67.07, "BMIM_C2":     15.04, "BMIM_C3":      57.12, "BMIM_COM":  0}
masses["BMIM_COM"] = masses["BMIM_C1"] +  masses["BMIM_C2"] +  masses["BMIM_C3"]
lj_sigmas       = {"PF6":       5.06, "BMIM_C1":        4.38, "BMIM_C2":      3.41, "BMIM_C3":       5.04, "BMIM_COM":  0}
lj_epsilons     = {"PF6":       2.56, "BMIM_C1":        2.56, "BMIM_C2":      0.36, "BMIM_C3":       1.83, "BMIM_COM":  0}
shortTypes  = ["A", "C1", "C2", "C3", "COM"]
lj_types       = ["PF6", "BMIM_C1", "BMIM_C2","BMIM_C3"]

cutoff_sigmafactor = 2.5
lj_cuts = {}
for t in lj_sigmas:
    lj_cuts[t] = cutoff_sigmafactor * lj_sigmas[t] 

#Num Particles and box
n_ionpairs= 300
n_part = n_ionpairs*2
density_si = 1.322 ;#g/cm^3, for bmimpf6, from paper; 1.322 (paper for 353K)
#density_si = 0.5 ;#g/cm^3, for bmimpf6, from paper; 1.322 (paper for 353K)
#rho_factor_bmim_pf6 = 0.003931
rho_factor_bmim_pf6 = 0.002
box_volume = n_ionpairs/rho_factor_bmim_pf6/density_si
box_l = box_volume**(1. / 3.)
print("Ion pairs:", n_ionpairs, "Box size:", box_l)
S.box_l=[box_l,box_l,box_l]
S.min_global_cut = 3.5

def combination_rule_epsilon(rule, eps1, eps2):
    if rule=="Lorentz":
        return (eps1*eps2)**0.5
    else:
        return ValueError("No combination rule defined")

def combination_rule_sigma(rule, sig1, sig2):
    if rule=="Berthelot":
        return (sig1+sig2)*0.5
    else:
        return ValueError("No combination rule defined")

# Lennard-Jones interactions parameters 
for i in range(len(lj_types)):
    for j in range(i, len(lj_types)):
        s=[lj_types[i], lj_types[j]]
        print("Mixing", s)
        lj_sig = combination_rule_sigma("Berthelot",lj_sigmas[s[0]], lj_sigmas[s[1]])
        lj_cut = combination_rule_sigma("Berthelot", lj_cuts[s[0]], lj_cuts[s[1]])
        lj_eps = combination_rule_epsilon("Lorentz", lj_epsilons[s[0]],lj_epsilons[s[1]])

        S.non_bonded_inter[types[s[0]], types[s[1]]].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

#Place Particles
rid=0
anion_ids=[]
cation_sites_ids=[]
cation_com_ids=[]
cation_c1_ids=[]
cation_c2_ids=[]
cation_c3_ids=[]

for i in range(n_ionpairs):
    S.part.add(id=rid, type=types["PF6"],  pos=np.random.random(3) * box_l, q=charges["PF6"], mass = masses["PF6"])
    anion_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1

for i in range(n_ionpairs):
    pos_com = np.random.random(3) * box_l
    S.part.add(id = rid, type = types["BMIM_COM"], pos = pos_com, mass = masses["BMIM_COM"], rinertia = [646.284, 585.158, 61.126], temp = 0, gamma = 0, rotation=[1,1,1])
    cation_com_ids.append(rid)
    com_id = rid
    rid+=1
    S.part.add(id = rid, type = types["BMIM_C1"], pos = pos_com + [0, -0.527, 1.365], q = charges["BMIM_C1"], virtual = 1)
    S.part[rid].vs_auto_relate_to(com_id)
    cation_c1_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1
    S.part.add(id = rid, type = types["BMIM_C2"], pos = pos_com + [0, 1.641, 2.987], q = charges["BMIM_C2"], virtual = 1)
    S.part[rid].vs_auto_relate_to(com_id)
    cation_c2_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1
    S.part.add(id = rid, type = types["BMIM_C3"], pos = pos_com + [0, 0.187, -2.389], q = charges["BMIM_C3"], virtual = 1)
    S.part[rid].vs_auto_relate_to(com_id)
    cation_c3_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1

#for i in cation_com_ids:
#    S.part[i].fix = [1,1,1]

#Energy minimization
Eminstart = time.time()
print("E before minimization:",S.analysis.energy()["total"])
#S.minimize_energy.init(f_max = 1000.0, gamma = 0.01, max_steps = 200, max_displacement= 0.01)
S.minimize_energy.init(f_max = 10.0, gamma = 0.01, max_steps = 10000, max_displacement= 0.01)
S.minimize_energy.minimize()
Emintime = (time.time()-Eminstart)
print("E after minimization:",S.analysis.energy()["total"], "in", Emintime, "s")
print("Min dist:", S.analysis.mindist([0,1],[0,1]))


#Forcecap equi
#print("Forcecap equi")
#cap = 10
#for i in xrange(100):
#    S.non_bonded_inter.set_force_cap(cap)
#    cap+=10
#    S.integrator.run(10)
#    S.galilei.kill_particle_motion()
#    print("E",S.analysis.energy()["total"])
#
#S.non_bonded_inter.set_force_cap(0)


#Actors
S.thermostat.set_langevin(kT=temperature_com, gamma=gamma_com)

#for i in xrange(1000):
#    S.integrator.run(10)


if args.gpup3m:
    print("Tune p3m GPU")
    p3m=P3M_GPU(bjerrum_length=l_b,accuracy=1e-4)
else:
    print("Tune p3m CPU")
    p3m=P3M(bjerrum_length=l_b,accuracy=1e-4)

S.actors.add(p3m)

if args.drude:
    print("Adding Drude bonds")
    #S.thermostat.turn_off()

    #Drude Bond
    drude_bond = DrudeBond(temp_com = temperature_com, gamma_com = gamma_com, temp_drude = temperature_drude, gamma_drude = gamma_drude, k = k_drude, r_cut = min(lj_sigmas.values())*0.5)
    S.bonded_inter.add(drude_bond)

    for i in anion_ids:
        add_drude_particle_to_core(S, S.part[i], drude_bond, i+1, types["PF6_D"], polarizations["PF6"], args.mass_drude, l_b)
    for i in cation_c1_ids:
        add_drude_particle_to_core(S, S.part[i], drude_bond, i+1, types["BMIM_C1_D"], polarizations["BMIM_C1"], args.mass_drude, l_b)
    for i in cation_c2_ids:
        add_drude_particle_to_core(S, S.part[i], drude_bond, i+1, types["BMIM_C2_D"], polarizations["BMIM_C2"], args.mass_drude, l_b)
    for i in cation_c3_ids:
        add_drude_particle_to_core(S, S.part[i], drude_bond, i+1, types["BMIM_C3_D"], polarizations["BMIM_C3"], args.mass_drude, l_b)

    if args.thole:
        add_all_thole(S)

    if args.intra_ex:
        #Setup bonds once
        setup_intramol_exclusion_bonds(S, [types["BMIM_C1_D"], types["BMIM_C2_D"], types["BMIM_C3_D"]], [types["BMIM_C1"], types["BMIM_C2"], types["BMIM_C3"]], [charges["BMIM_C1"], charges["BMIM_C2"], charges["BMIM_C3"]])

        #add sr ex bonds per molecule
        for i in cation_c1_ids:
            add_intramol_exclusion_bonds(S, [i+1,i+3,i+5], [i,i+2,i+4])



#icc=ICC(n_icc=nicc_tot , convergence=1e-6, relaxation=0.75, ext_field=[0,0,0], max_iterations=100, first_id=0, eps_out=1, normals=iccNormals, areas=iccAreas, sigmas=iccSigmas, epsilons=iccEpsilons)
#S.actors.add(icc)


if args.timing:
    print("Timing")
    start = time.time()
    S.integrator.run(1000)
    time_per_step = (time.time()-start)/1000.0
    ns_per_hour = 3600.0 * time_step_fs * 1e-6 / time_per_step 
    ns_per_day = 24.0 * ns_per_hour
    print("Yield:", ns_per_day, "ns/day")
    exit()

if args.visu:
    def main():
        print("\n--->Integration")

        while (True):
            S.integrator.run(10)
            visualizer.update()

    #Start simulation in seperate thread
    t = Thread(target=main)
    t.daemon = True
    t.start()

    #Start blocking visualizer
    visualizer.start()

else:
    print("\nEquilibration\n")
    n_int_steps = 50
    n_int_cycles = 9000

    for i in range(n_int_cycles):
        S.integrator.run(n_int_steps)
        print(100.0*i/n_int_cycles, "% Done")

    print("\nIntegrate\n")

    obs_Temps = []
    obs_Ekin = []

    obs_Etot = []
    obs_Momentum = []

    obs_dists_drude_vec = []

    n_int_steps = 100
    n_int_cycles = int(args.walltime * 3600.0 / time_per_step / n_int_steps)
    print("Performing", n_int_cycles, "x", n_int_steps, "steps, which is", args.walltime*ns_per_hour, "ns simulation time")

    file_traj = open(path + "traj.xyz", "w")

    for i in range(n_int_cycles):
        S.integrator.run(n_int_steps)
      
        file_traj.write(str(n_part*2) + '\n')
        file_traj.write("t(ns) = " + str(time_step_ns*i*n_int_steps)+'\n')
        for j in range(n_part):
            if args.drude:
                #Drude elongation
                vec = S.part[2*j].pos - S.part[2*j+1].pos
                obs_dists_drude_vec.append([vec[0],vec[1],vec[2], np.linalg.norm(vec)])
            #Core trajectories
                file_traj.write(shortTypes[S.part[2*j+1].type] + " " + ' '.join(map(str,S.part[2*j+1].pos_folded)) + '\n')
            file_traj.write(shortTypes[S.part[id_f*j].type] + " " + ' '.join(map(str,S.part[id_f*j].pos_folded)) + '\n')

        #Ana append for rdfs
        S.analysis.append()

        if i%100==0:
            print(100*i/n_int_cycles, "% Done   ")
        

    #Temp/Ekin histos
    if args.drude:
        Ekin_com_histo = np.histogram(np.array(obs_Ekin)[:,0], bins = 100)
        temp_com_histo = np.histogram(np.array(obs_Temps)[:,0], bins = 100)
        np.savetxt(path + "Ekin_com.dat", np.vstack((Ekin_com_histo[1][:-1],Ekin_com_histo[0])).T)
        np.savetxt(path + "temp_com.dat", np.vstack((temp_com_histo[1][:-1],temp_com_histo[0])).T)

        Ekin_dist_histo = np.histogram(np.array(obs_Ekin)[:,2], bins = 100)
        temp_dist_histo = np.histogram(np.array(obs_Temps)[:,2], bins = 100)
        np.savetxt(path + "Ekin_dist.dat", np.vstack((Ekin_dist_histo[1][:-1],Ekin_dist_histo[0])).T)
        np.savetxt(path + "temp_dist.dat", np.vstack((temp_dist_histo[1][:-1],temp_dist_histo[0])).T)

    Ekin_core_histo = np.histogram(np.array(obs_Ekin)[:,1], bins = 100)
    temp_core_histo = np.histogram(np.array(obs_Temps)[:,1], bins = 100)
    np.savetxt(path + "Ekin_core.dat", np.vstack((Ekin_core_histo[1][:-1],Ekin_core_histo[0])).T)
    np.savetxt(path + "temp_core.dat", np.vstack((temp_core_histo[1][:-1],temp_core_histo[0])).T)

    #Global
    np.savetxt(path + "energies.dat", obs_Etot)
    np.savetxt(path + "momentum.dat", obs_Momentum)

    if args.drude:
        #Drude histos
        drude_dist_x_histo = np.histogram(np.array(obs_dists_drude_vec)[:,0], bins = 100)
        drude_dist_y_histo = np.histogram(np.array(obs_dists_drude_vec)[:,1], bins = 100)
        drude_dist_z_histo = np.histogram(np.array(obs_dists_drude_vec)[:,2], bins = 100)
        drude_dist_len_histo = np.histogram(np.array(obs_dists_drude_vec)[:,3], bins = 100)

        np.savetxt(path + "drude_dist_x.dat", np.vstack((drude_dist_x_histo[1][:-1],drude_dist_x_histo[0])).T)
        np.savetxt(path + "drude_dist_y.dat", np.vstack((drude_dist_y_histo[1][:-1],drude_dist_y_histo[0])).T)
        np.savetxt(path + "drude_dist_z.dat", np.vstack((drude_dist_z_histo[1][:-1],drude_dist_z_histo[0])).T)
        np.savetxt(path + "drude_dist_len.dat", np.vstack((drude_dist_len_histo[1][:-1],drude_dist_len_histo[0])).T)

    #RDFs
    rdf_bins = 100
    r_min  = 0.0
    r_max  = S.box_l[0]/2.0
    r,rdf_00 = S.analysis.rdf(rdf_type='<rdf>', 
                                type_list_a=[types["Anion"]],
                                type_list_b=[types["Anion"]], 
                                r_min=r_min,
                                r_max=r_max, 
                                r_bins=rdf_bins)

    r,rdf_11 = S.analysis.rdf(rdf_type='<rdf>', 
                                type_list_a=[types["Cation"]],
                                type_list_b=[types["Cation"]], 
                                r_min=r_min,
                                r_max=r_max, 
                                r_bins=rdf_bins)
    r,rdf_01 = S.analysis.rdf(rdf_type='<rdf>', 
                                type_list_a=[types["Anion"]],
                                type_list_b=[types["Cation"]], 
                                r_min=r_min,
                                r_max=r_max, 
                                r_bins=rdf_bins)

    rdf_fp = open(path + "rdf.dat", 'w')
    for i in range(rdf_bins):
        rdf_fp.write("%1.5e %1.5e %1.5e %1.5e\n" % (r[i], rdf_00[i], rdf_11[i], rdf_01[i]))
    rdf_fp.close()

    #MSD
    file_traj.close()

    print("Done. Output path:",path)

    #msd=[]
    #for t in range(len(cores_t)):
    #    msd.append([t*time_step_ns*n_int_steps, ((np.array(cores_t[t])-np.array(cores_t[0]))**2).sum()/n_part])

    #np.savetxt(path + "msd.dat", msd)

