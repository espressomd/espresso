from __future__ import print_function
import espressomd
from espressomd.visualization_opengl import openGLLive
import numpy as np
from math import *

system = espressomd.System()
box_l = 15
system.box_l = [box_l, box_l, box_l]

# cup
cup_top_circ = 21
cup_bot_circ = 15
cup_height = 6
for i in range(cup_height):
    circ = cup_bot_circ + i*(cup_top_circ-cup_bot_circ)/float(cup_height-1)
    rad = circ/(2.0*np.pi)
    alpha = 2.0*np.pi/int(circ)
    posy = 2.0+i
    for j in range(int(circ)):
        posx = box_l/2.0 + rad * sin(j*alpha+(np.pi/2.0))
        posz = box_l/2.0 + rad * cos(j*alpha+(np.pi/2.0))
        system.part.add(pos = [posx, posy, posz], type = 0)

# cup bottom
rad  = cup_bot_circ/(2.0*np.pi)
posy = 2.0
while (rad > 1.0):
    rad -= 0.9
    circ = 2.0*np.pi*rad
    alpha = 2.0*np.pi/int(circ)
    for j in range(int(circ)):
        posx = box_l/2.0 + rad*sin(j*alpha+(np.pi/2.0))
        posz = box_l/2.0 + rad*cos(j*alpha+(np.pi/2.0))
        system.part.add(pos = [posx, posy, posz], type = 0)


# cup handle
hand_rad = (cup_height-4.0)/sqrt(2.0)
hand_circ = (1.5*np.pi*hand_rad)
hand_xoff = (cup_bot_circ+cup_top_circ)/(4.0*np.pi)+1.2
hand_yoff = 2.0+cup_height/2.0 -0.2
alpha = 2.0*np.pi/int(4.0*hand_circ/3.0)
beta =  sin((cup_top_circ-cup_bot_circ)/(2.0*np.pi*cup_height-1))
beta = beta - np.pi/8.0
posz = (box_l/2.0)+0.5
for i in range(int(hand_circ)):
    posx = hand_xoff + box_l/2.0 + hand_rad*sin(i*alpha+beta)
    posy = hand_yoff + hand_rad*cos(i*alpha+beta)
    system.part.add(pos = [posx, posy, posz], type = 0)

# saucer
saucer_circ = 30
s_rad_o = saucer_circ/(2.0*np.pi)
s_rad_i = cup_bot_circ/(2.0*np.pi)
n_saucer = int(s_rad_o-s_rad_i) + 1
n_ci = 0
for i in range(n_saucer):
    n_ci += int(saucer_circ-(i*2.0*np.pi))

ci_val = -len(system.part)/float(n_ci)
for i in range(n_saucer):
    rad = s_rad_o-i
    alpha = 2.0*np.pi/int(saucer_circ-(i*2.0*np.pi))
    posy = 2.3 - 0.5*i
    for j in range(int(saucer_circ-(i*2.0*np.pi))):
        posx = box_l/2.0 + rad*sin(j*alpha)
        posz =  box_l/2.0 + rad*cos(j*alpha)
        system.part.add(pos = [posx, posy, posz], type = 1)

# steam
fene = espressomd.interactions.FeneBond(k=10.0, d_r_max=2.0, r_0=0.1)
system.bonded_inter.add(fene)

n_steam = 10
l_steam = 10
rad = (cup_top_circ-12.0)/(2.0*np.pi)
alpha = 2.0*np.pi/int(n_steam)
for i in range(n_steam):
    for j in range(l_steam):
        posx = box_l/2.0 + rad*sin(i*alpha + j*0.6)
        posz = box_l/2.0 + rad*cos(i*alpha + j*0.6)
        posy = cup_height + j * 0.1 *rad
        system.part.add(pos = [posx, posy, posz], type = 2)
        pid = len(system.part)-1

        if j == 0:
            system.part[pid].fix = [1,1,1]
        else:
            system.part[pid].bonds = (fene, pid-1)

        if j == l_steam-1:
            system.part[pid].ext_force = [0,5.0,0]


system.time_step = 0.00022
system.cell_system.skin = 0.4

system.thermostat.set_langevin(kT=0.0, gamma=0.07)
WCA_cut = 2.**(1. / 6.)

lj_eps = 1.0
lj_sig = 0.7
lj_cut = WCA_cut * lj_sig
for i in range(2):
    for j in range(i,2):
        system.non_bonded_inter[i,j].lennard_jones.set_params(epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

visualizer = openGLLive(system, 
						background_color = [0.2,0.2,0.3],
                        particle_sizes = [0.6,0.7,0.3], 
						particle_type_materials = ['bright', 'gold', 'chrome'],
						bond_type_materials = ['chrome'],
                        particle_type_colors = [[0, 0, 0.8, 1], [0.8, 0, 0, 1], [0.8, 0.8, 0.8, 1]],
                        bond_type_colors = [[0.2, 0.2, 0.2, 1], [0, 1, 1, 1], [1, 0.5, 0, 1]],
                        bond_type_radius = [0.1],
                        spotlight_brightness = 1.0,
                        light_brightness = 1.0,
                        ext_force_arrows = False,
                        draw_axis = False,
                        draw_box = False,
                        drag_enabled = True )

def rotate():
	visualizer.camera.rotateSystemXL()

visualizer.registerCallback(rotate, interval = 16)

visualizer.run(1)
