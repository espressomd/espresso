#!/usr/bin/python

from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
from espressomd import minimize_energy
from espressomd.interactions import *
import numpy
from threading import Thread
from math import *
from espressomd.visualization_opengl import *
from espressomd.shapes import *

print('8Ball BILLARD - An Espresso Visualizer Demo\nControls:\nNumpad 4/6: Adjust Angle\nNumpad 2/8: Adjust Impulse\nNumpad 5: Shoot')

#ESPRESSO
system = espressomd.System()
table_dim = [2.24,1.12]
system.box_l = [table_dim[0], 3, table_dim[1]]

visualizer = openGLLive(system,
        ext_force_arrows = True, 
        ext_force_arrows_scale =  [0.02], 
        background_color = [0.5,0.4,0.5], 
        drag_enabled = False, 
        constraint_type_colors = [[0.039,0.424,0.011,1.0]],
        particle_type_colors = [[1,1,1,1],[1,0,0,1],[0,0,1,1],[0.2,0.2,0.2,1]],
        camera_position = [ 1.12, 2.8, 0.56],
        window_size = [1000,600],
        draw_axis = False,
        light_brightness = 5.0)

stopped = True
angle = numpy.pi*0.5
impulse = 10.0

def decreaseAngle():
    global angle,impulse
    if stopped:
        angle += 0.01
        system.part[0].ext_force = impulse*np.array([sin(angle),0,cos(angle)])

def increaseAngle():
    global angle,impulse
    if stopped:
        angle -= 0.01
        system.part[0].ext_force = impulse*np.array([sin(angle),0,cos(angle)])

def decreaseImpulse():
    global impulse,angle
    if stopped:
        impulse -= 0.5
        system.part[0].ext_force = impulse*np.array([sin(angle),0,cos(angle)])

def increaseImpulse():
    global impulse,angle
    if stopped:
        impulse += 0.5
        system.part[0].ext_force = impulse*np.array([sin(angle),0,cos(angle)])

def fire():
    global stopped
    if stopped:
        stopped = False
        system.part[0].v = system.part[0].v + impulse * np.array([sin(angle),0,cos(angle)])
        system.part[0].fix = [0,1,0]
        system.part[0].ext_force = [0,0,0]

visualizer.keyboardManager.registerButton(KeyboardButtonEvent('4',KeyboardFireEvent.Hold,decreaseAngle))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('6',KeyboardFireEvent.Hold,increaseAngle))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('2',KeyboardFireEvent.Hold,decreaseImpulse))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('8',KeyboardFireEvent.Hold,increaseImpulse))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('5',KeyboardFireEvent.Pressed,fire))

def main():
    global stopped

    system.time_step = 0.00008
    system.cell_system.skin = 0.4

    table_h = 0.5
    ball_diam = 0.0572
    hole_dist = 0.02
    hole_rad = 0.08
    hole_score_rad = 0.1
    hole_pos = [[hole_dist, table_h, hole_dist],
                [hole_dist, table_h, table_dim[1] - hole_dist],
                [table_dim[0] - hole_dist, table_h, hole_dist],
                [table_dim[0] - hole_dist, table_h, table_dim[1] - hole_dist],
                [table_dim[0]*0.5, table_h, table_dim[1] - hole_dist],
                [table_dim[0]*0.5, table_h, hole_dist]]
    types = {'cue_ball': 0,'striped_ball':1, 'solid_ball':2,'black_ball':3, 'table':4,'wall':5,'hole':6}


    system.constraints.add(shape=Wall(dist=table_h,normal=[0.0,1.0,0.0]),particle_type=types['table'],penetrable=1)

    system.constraints.add(shape=Wall(dist=0.01,normal=[1.0,0.0,0.0]),particle_type=types['wall'],penetrable=1)
    system.constraints.add(shape=Wall(dist=-(table_dim[0]-0.01),normal=[-1.0,0.0,0.0]),particle_type=types['wall'],penetrable=1)
    system.constraints.add(shape=Wall(dist=0.01,normal=[0.0,0.0,1.0]),particle_type=types['wall'],penetrable=1)
    system.constraints.add(shape=Wall(dist=-(table_dim[1]-0.01),normal=[0.0,0.0,-1.0]),particle_type=types['wall'],penetrable=1)
    for h in hole_pos:
        system.constraints.add(shape=Cylinder(center=(np.array(h)-np.array([0,table_h*0.5,0])).tolist(), axis=[0,1,0],radius = hole_rad, length = 1.02*table_h, direction = 1),particle_type=types['hole'], penetrable=1)

    lj_eps = np.array([1])
    lj_sig = np.array([ball_diam])
    lj_cut = lj_sig*2.0**(1.0/6.0)
    lj_cap = 20
    mass = np.array([0.17])

    num_types=len(lj_sig)

    #LENNARD JONES
    def mix_eps(eps1,eps2,rule='LB'):
        return sqrt(eps1*eps2)
    def mix_sig(sig1,sig2,rule='LB'):
        return 0.5*(sig1+sig2)
    for t1 in range(4):
        for t2 in range(6):
            system.non_bonded_inter[t1, t2].lennard_jones.set_params(
                epsilon=mix_eps(lj_eps[0],lj_eps[0]), 
                sigma=mix_sig(lj_sig[0],lj_sig[0]),
                cutoff=mix_sig(lj_cut[0],lj_cut[0]), 
                shift="auto")

    ball_y = table_h+ball_diam*1.5

    #PARTICLES
    ball_start_pos = [table_dim[0]*0.25, ball_y, table_dim[1]*0.5]
    system.part.add(id=0, pos=ball_start_pos ,type=types['cue_ball'],mass=mass[0])
    spawnpos = []
    spawnpos.append(ball_start_pos)
    ball = system.part[0]

    d = lj_sig[0]*1.15
    a1 = np.array([d*sqrt(3)/2.0,0, -0.5*d])
    a2 = np.array([d*sqrt(3)/2.0,0, 0.5*d])
    sp = [system.box_l[0]*0.7, ball_y, system.box_l[2]*0.5+lj_sig[0]*0.5]
    pid = 1
    order = [
    types['solid_ball'],
    types['striped_ball'],types['solid_ball'],
    types['solid_ball'],types['black_ball'],types['striped_ball'],
    types['striped_ball'], types['solid_ball'],types['striped_ball'], types['solid_ball'],
    types['solid_ball'], types['striped_ball'], types['striped_ball'],types['solid_ball'], types['striped_ball']]

    for i in range(5):
        for j in range(i+1):
            N=i+1
            t = order[pid-1]
            pos = sp + a1*(N-j) + a2*j
            system.part.add(id = pid, pos=pos, mass=mass[0],type =t, fix = [0,1,0])
            spawnpos.append(pos)
            pid += 1

    ball.ext_force = impulse*np.array([sin(angle),0,cos(angle)])
    ball.fix = [1,1,1]
    system.thermostat.set_langevin(kT=0, gamma=0.8)
    #ELECTROSTATICS
    #	p3m = electrostatics.P3M(prefactor=50, accuracy=1e-2)
    #	system.actors.add(p3m)

    cleared_balls= [0,0]
    while True:	
        system.integrator.run(1)

        vsum = 0
        for p in system.part:
            vsum += numpy.linalg.norm(p.v)

            for h in hole_pos:
                
                d=((p.pos_folded[0]-h[0])**2 + (p.pos_folded[2]-h[2])**2)**0.5
                if (d < hole_score_rad):
                    if p.id == 0:
                        p.pos = ball_start_pos
                        p.v = [0,0,0]
                    elif p.id == 5:
                        for p in system.part:
                            p.pos = spawnpos[p.id]
                            p.v = [0,0,0]
                            p.fix = [0,1,0]
                        ball.fix = [1,1,1]
                        ball.ext_force = impulse*np.array([sin(angle),0,cos(angle)])
                        stoppen = True
                    else:
                        t = p.type-1 
                        cleared_balls[t] += 1
                        if t == 0:
                            z=table_dim[1]-lj_sig[0]*0.6
                        else:
                            z=lj_sig[0]*0.6
                        p.pos = [cleared_balls[t] * lj_sig[0] * 1.5 , 1.1, z]
                        p.fix = [1,1,1]
                        p.v = [0,0,0]

        if not stopped and vsum < 0.3:
            stopped = True
            ball.fix = [1,1,1]
            for p in system.part:
                p.v = [0,0,0]
            ball.ext_force = impulse*np.array([sin(angle),0,cos(angle)])

        visualizer.update()


t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
