# Copyright (C) 2010-2019 The ESPResSo project
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
"""
Game based on Maxwell's demon, a thought experiment used to teach statistical
thermodynamics. The user has to scoop particles from a chamber and guide them
to another chamber through a channel with the help of a snake controlled by a
gamepad or the keyboard. The particle imbalance between chambers creates
a pressure gradient that makes it harder to move particles to the chamber
with an excess of particles.
"""

from threading import Thread
import numpy as np
import time

import espressomd
import espressomd.shapes
import espressomd.minimize_energy
from espressomd.visualization_opengl import openGLLive, KeyboardButtonEvent, KeyboardFireEvent

required_features = ["LENNARD_JONES", "WCA", "MASS",
                     "EXTERNAL_FORCES", "LANGEVIN_PER_PARTICLE"]
espressomd.assert_features(required_features)

print("""THE CHAMBER GAME

YOUR GOAL IS TO SCOOP ALL BLUE PARTICLES INTO THE RIGHT BOX.
GREEN/RED SPHERES CAN BE PICKED UP TO INCREASE/DECREASE
THE TEMPERATURE IN THE CHAMBER WHERE THEY ARE COLLECTED.""")

try:
    import pygame
    has_pygame = True
    print("\nCONTROLS:"
          "\nMOVE: (JOYSTICK AXIS), (KEYBOARD i/j/k/l)"
          "\nACTION BUTTON: (JOYSTICK A), (KEYBOARD p)"
          "\nRESTART: (JOYSTICK START), (KEYBOARD b)")
except BaseException:
    has_pygame = False
    print("\nCONTROLS:"
          "\nMOVE: (KEYBOARD i/j/k/l)"
          "\nACTION BUTTON: (KEYBOARD p)"
          "\nRESTART: (KEYBOARD b)")

box = np.array([1500.0, 500.0, 150.0])
system = espressomd.System(box_l=box)

# PARAMETERS

# PHYSICS
temperature_snake = 0.0
gamma_snake_head = 1.0
gamma_snake_bead = 15.0

temperature_bubbles = 10000.0
temp_l = temperature_bubbles
temp_r = temperature_bubbles
temp_max = 1e5
gamma_bubbles = 0.5

temperature = 1.0
gamma = 1.0
system.time_step = 0.001

# SNAKE
snake_n = 10
snake_head_sigma = 50.0
snake_bead_sigma = 20.0
snake_length = (snake_n - 1) * snake_bead_sigma + snake_head_sigma
snake_startpos = [snake_head_sigma, box[1] - snake_head_sigma, box[2] * 0.5]
snake_head_type = 0
snake_bead_type = 1
snake_head_mass = 1000.0
snake_bead_mass = 10.0
harmonic_k = 500.0 * snake_bead_mass

# PORE
pore_length = box[0] * 0.25
pore_xl = box[0] * 0.5 - pore_length * 0.5
pore_xr = box[0] * 0.5 + pore_length * 0.5
cylinder_type = 2
cylinder_sigma = 1.0
pore_radius = snake_head_sigma * 1.3

# CONTROL
move_force = 70000.0
expl_range = 200.0
expl_force = 20000.0

# BUBBLES
bubble_type = 3
bubble_sigma = 36.0
bubble_snake_eps = 10
bubble_bubble_eps = 10000.0
bubble_mass = 50.0
bubbles_n = 180

# TEMP CHANGE PARTICLE
temp_change_radius = 25
temp_change_inc_type = 4
temp_change_dec_type = 5
dtemp = 1000.0

# VISUALIZER
zoom = 10

visualizer = openGLLive(
    system,
    window_size=[800, 600],
    draw_axis=False,
    particle_sizes=[
        snake_head_sigma * 0.5,
        snake_bead_sigma * 0.5,
        cylinder_sigma,
        bubble_sigma * 0.5,
        temp_change_radius,
        temp_change_radius],
    particle_type_colors=[[1, 1, 0],
                          [1, 0, 1],
                          [0, 0, 1],
                          [0, 1, 1],
                          [0, 1, 0],
                          [1, 0, 0],
                          [0.5, 0, 1]],
    constraint_type_colors=[[1, 1, 1]],
    camera_position=[snake_startpos[0],
                     snake_startpos[1],
                     system.box_l[2] * zoom],
    camera_target=snake_startpos)


# JOYPAD CONTROL
if has_pygame:
    pygame.init()
    pygame.joystick.init()

    # CHECK FOR JOYSTICKS
    if pygame.joystick.get_count() > 0:
        joystick = pygame.joystick.Joystick(0)
        joystick.init()
        joystick_control = True
    else:
        joystick_control = False

# CELLSYSTEM

system.cell_system.skin = 3.0
system.cell_system.set_domain_decomposition(use_verlet_lists=False)

# BONDS

harmonic_head = espressomd.interactions.HarmonicBond(
    k=harmonic_k, r_0=0.5 * (snake_head_sigma + snake_bead_sigma))
harmonic_bead = espressomd.interactions.HarmonicBond(
    k=harmonic_k, r_0=snake_bead_sigma)
system.bonded_inter.add(harmonic_head)
system.bonded_inter.add(harmonic_bead)

# PARTICLES

# SNAKE

for i in range(snake_n):
    if i == 0:
        p_head = system.part.add(
            pos=snake_startpos,
            type=snake_head_type,
            fix=[False, False, True],
            mass=snake_head_mass,
            temp=temperature_snake,
            gamma=gamma_snake_head)
    else:
        system.part.add(
            pos=snake_startpos
            + np.array([0, -1, 0])
            * (0.5 * (snake_head_sigma + snake_bead_sigma)
               + (i - 1) * snake_bead_sigma),
            bonds=(harmonic_bead if (i > 1) else harmonic_head, i - 1),
            type=snake_bead_type,
            fix=[False, False, True],
            mass=snake_bead_mass,
            temp=temperature_snake,
            gamma=gamma_snake_bead)

# NB INTER

WCA_cut = 2.0**(1. / 6.)

system.non_bonded_inter[snake_head_type, snake_head_type].wca.set_params(
    epsilon=1.0, sigma=snake_head_sigma)

sm = 0.5 * (snake_head_sigma + snake_bead_sigma)
system.non_bonded_inter[snake_bead_type, snake_head_type].wca.set_params(
    epsilon=1.0, sigma=sm)

system.non_bonded_inter[snake_bead_type, snake_bead_type].wca.set_params(
    epsilon=1.0, sigma=snake_bead_sigma)

sm = 0.5 * (snake_head_sigma + cylinder_sigma)
system.non_bonded_inter[snake_head_type, cylinder_type].wca.set_params(
    epsilon=10.0, sigma=sm)

sm = 0.5 * (snake_bead_sigma + cylinder_sigma)
system.non_bonded_inter[snake_bead_type, cylinder_type].wca.set_params(
    epsilon=10.0, sigma=sm)

sm = 0.5 * (bubble_sigma + snake_bead_sigma)
system.non_bonded_inter[snake_bead_type, bubble_type].wca.set_params(
    epsilon=bubble_snake_eps, sigma=sm)

sm = 0.5 * (bubble_sigma + snake_head_sigma)
system.non_bonded_inter[snake_head_type, bubble_type].wca.set_params(
    epsilon=1.0, sigma=sm)

sm = 0.5 * (bubble_sigma + cylinder_sigma)
system.non_bonded_inter[bubble_type, cylinder_type].lennard_jones.set_params(
    epsilon=1000.0, sigma=sm, cutoff=2.5 * sm, shift="auto")

system.non_bonded_inter[bubble_type, bubble_type].lennard_jones.set_params(
    epsilon=bubble_bubble_eps, sigma=bubble_sigma, cutoff=2.5 * bubble_sigma, shift="auto")

# CONSTRAINTS

system.constraints.add(shape=espressomd.shapes.Wall(
    dist=0, normal=[1, 0, 0]), particle_type=cylinder_type, penetrable=True)
system.constraints.add(shape=espressomd.shapes.Wall(
    dist=-box[0], normal=[-1, 0, 0]), particle_type=cylinder_type, penetrable=True)
system.constraints.add(shape=espressomd.shapes.Wall(
    dist=0, normal=[0, 1, 0]), particle_type=cylinder_type, penetrable=True)
system.constraints.add(shape=espressomd.shapes.Wall(
    dist=-box[1], normal=[0, -1, 0]), particle_type=cylinder_type, penetrable=True)

system.constraints.add(shape=espressomd.shapes.SimplePore(
    center=0.5 * box, axis=[1, 0, 0], length=pore_length, radius=pore_radius,
    smoothing_radius=5), particle_type=cylinder_type, penetrable=True)


# BUBBLES
n = 0

while n < bubbles_n:
    # bpos = [pore_xr +  np.random.random() * (pore_xr - pore_xl -
    # snake_head_sigma*4) + snake_head_sigma * 2, np.random.random() * box[1],
    # box[2]*0.5]
    bpos = [np.random.random() * (pore_xl - snake_head_sigma * 4) +
            snake_head_sigma * 2, np.random.random() * box[1], box[2] * 0.5]
    system.part.add(
        pos=bpos,
        type=bubble_type,
        fix=[False, False, True],
        mass=bubble_mass,
        temp=temperature_bubbles,
        gamma=gamma_bubbles)
    testid = len(system.part) - 1
    n += 1

    if np.min([system.distance(system.part[testid], p.pos)
               for p in system.part if p.id != testid]) < bubble_sigma * 0.5:
        system.part[testid].remove()
        n -= 1

p_bubbles = np.where(system.part[:].type == bubble_type)[0]

# TEMP CHANGE PARTICLES
bpos = [np.random.random() * (pore_xl - snake_head_sigma * 4) +
        snake_head_sigma * 2, np.random.random() * box[1], box[2] * 0.5]
p_temp_inc = system.part.add(
    pos=bpos,
    type=temp_change_inc_type,
    fix=[True, True, True])

bpos = [pore_xr
        + np.random.random() * (pore_xr - pore_xl - snake_head_sigma * 4)
        + snake_head_sigma * 2,
        np.random.random() * box[1],
        box[2] * 0.5]
p_temp_dec = system.part.add(
    pos=bpos,
    type=temp_change_dec_type,
    fix=[True, True, True])

# MINIMIZE ENERGY

energy = system.analysis.energy()
#print("Before Minimization: E_total = {}".format(energy['total']))
espressomd.minimize_energy.steepest_descent(system, f_max=100, gamma=30.0,
                                            max_steps=10000,
                                            max_displacement=0.01)
energy = system.analysis.energy()
#print("After Minimization: E_total = {}".format(energy['total']))

p_startpos = system.part[:].pos

# THERMOSTAT
system.thermostat.set_langevin(kT=temperature, gamma=gamma, seed=42)

# CONTROL CALLBACKS
F_act_k = np.zeros(2)
F_act_j = np.zeros(2)


def move_up_set():
    global F_act_k
    F_act_k[1] = 1.0
    set_particle_force()


def move_down_set():
    global F_act_k
    F_act_k[1] = -1.0
    set_particle_force()


def move_updown_reset():
    global F_act_k
    F_act_k[1] = 0
    set_particle_force()


def move_left_set():
    global F_act_k
    F_act_k[0] = -1.0
    set_particle_force()


def move_right_set():
    global F_act_k
    F_act_k[0] = 1.0
    set_particle_force()


def move_leftright_reset():
    global F_act_k
    F_act_k[0] = 0
    set_particle_force()


def set_particle_force():
    global F_act_j, F_act_k
    F_control_tot = np.append(np.clip(F_act_k + F_act_j, -1, 1), 0)
    system.part[0].ext_force = move_force * F_control_tot


def restart():
    system.part[:].pos = p_startpos
    system.galilei.kill_particle_motion()
    system.galilei.kill_particle_forces()


expl_time = 0
exploding = False


def explode():
    global exploding, expl_time
    if not exploding:
        exploding = True
        expl_time = time.time()
        for p in system.part[p_bubbles]:
            dv = p.pos - p_head.pos
            lv = np.linalg.norm(dv)
            if lv < expl_range:
                p.v = dv / lv / lv * expl_force


# KEYBOARD CONTROLS
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('i', KeyboardFireEvent.Pressed, move_up_set))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('k', KeyboardFireEvent.Pressed, move_down_set))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('i', KeyboardFireEvent.Released, move_updown_reset))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('k', KeyboardFireEvent.Released, move_updown_reset))

visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('j', KeyboardFireEvent.Pressed, move_left_set))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('l', KeyboardFireEvent.Pressed, move_right_set))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('j', KeyboardFireEvent.Released, move_leftright_reset))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('l', KeyboardFireEvent.Released, move_leftright_reset))

visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('p', KeyboardFireEvent.Pressed, explode))

visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('b', KeyboardFireEvent.Pressed, restart))

# MAIN LOOP


def main():
    global F_act_j, F_act_k, temp_l, temp_r, exploding, expl_time

    def T_to_g(temp):
        return 0.1 + 5.0 / (1.0 + 0.001 * temp)

    zoom_eq = 5.0
    zoom_v = 0.0
    zoom_a = 0.0
    zoom = zoom_eq
    zoom_dt = 0.01

    ud_cnt = 0
    tincF = 0
    tdecF = 0
    exploding = False
    button_A_old = 0
    button_Start_old = 0
    while True:

        # INTEGRATE
        system.integrator.run(1)

        if p_head.pos[0] > pore_xl and p_head.pos[0] < pore_xr:
            z_eq = 10.0
            v_f = 0.1
        else:
            z_eq = zoom_eq
            v_f = 1.0

        # CAMERA TRACKING
        zoom_a = (z_eq - zoom) * 0.2 - zoom_v * 0.8 + v_f * \
            0.005 * np.linalg.norm(system.part[0].v)
        zoom_v += zoom_a * zoom_dt
        zoom += zoom_v * zoom_dt + zoom_a * zoom_dt * zoom_dt
        camPos = np.copy(system.part[0].pos) - box * 0.5
        camPos[2] = box[2] * zoom
        camTarget = system.part[0].pos - box * 0.5
        t = camPos - camTarget
        r = np.linalg.norm(t)
        visualizer.camera.state_pos = camPos
        visualizer.camera.state_target = -t / r
        visualizer.camera.update_modelview()

        # COUNT L/R
        ud_cnt += 1
        if ud_cnt > 100:
            ud_cnt = 0
            pl = system.part.select(
                lambda p: p.pos[0] < pore_xl and p.type == bubble_type)
            pr = system.part.select(
                lambda p: p.pos[0] > pore_xr and p.type == bubble_type)
            Nl = len(pl)
            Nr = len(pr)
            for p in pl:
                p.temp = temp_l
                p.gamma = T_to_g(temp_l)
            for p in pr:
                p.temp = temp_r
                p.gamma = T_to_g(temp_r)

            w = visualizer.specs['window_size']
            visualizer.user_texts = [
                [[20, w[1] - 20], 'LEFT: {}   RIGHT: {}'.format(Nl, Nr)],
                [[20, w[1] - 40], 'TEMPERATURE LEFT: {:.0f}   TEMPERATURE RIGHT: {:.0f}'.format(temp_l, temp_r)]]
            # [[w[0] * 0.5, w[1] - 60], 'GAMMA LEFT: {:0.4f}   GAMMA RIGHT: {:0.4f}'.format( T_to_g(temp_l), T_to_g(temp_r))]]

        # TEMP CHANGE COLLISION
        repos_temp_inc = False
        repos_temp_dec = False
        if np.linalg.norm(
                p_head.pos - p_temp_inc.pos) < temp_change_radius + snake_head_sigma * 0.5:
            repos_temp_inc = True
            if p_temp_inc.pos[0] > box[0] * 0.5:
                temp_r += dtemp
                if temp_r > temp_max:
                    temp_r = temp_max
            else:
                temp_l += dtemp
                if temp_l > temp_max:
                    temp_l = temp_max
        if np.linalg.norm(
                p_head.pos - p_temp_dec.pos) < temp_change_radius + snake_head_sigma * 0.5:
            repos_temp_dec = True
            if p_temp_dec.pos[0] > box[0] * 0.5:
                temp_r -= dtemp
                if temp_r < 0:
                    temp_r = 0.0
                for p in system.part[p_bubbles]:
                    if p.pos[0] > pore_xr:
                        p.v = [0, 0, 0]
            else:
                temp_l -= dtemp
                if temp_l < 0:
                    temp_l = 0.0
                for p in system.part[p_bubbles]:
                    if p.pos[0] < pore_xl:
                        p.v = [0, 0, 0]

        # PLACE TEMP CHANGE PARTICLES
        tincF += 1
        tdecF += 1
        if repos_temp_inc or tincF > 5000:
            tincF = 0
            if np.random.random() < 0.5:
                p_temp_inc.pos = [np.random.random()
                                  * (pore_xl - snake_head_sigma * 4)
                                  + snake_head_sigma * 2,
                                  np.random.random() * box[1],
                                  box[2] * 0.5]
            else:
                p_temp_inc.pos = [pore_xr
                                  + np.random.random()
                                  * (pore_xr - pore_xl - snake_head_sigma * 4)
                                  + snake_head_sigma * 2,
                                  np.random.random() * box[1],
                                  box[2] * 0.5]

        if repos_temp_dec or tdecF > 5000:
            tdecF = 0
            if np.random.random() < 0.5:
                p_temp_dec.pos = [np.random.random()
                                  * (pore_xl - snake_head_sigma * 4)
                                  + snake_head_sigma * 2,
                                  np.random.random() * box[1],
                                  box[2] * 0.5]
            else:
                p_temp_dec.pos = [pore_xr
                                  + np.random.random()
                                  * (pore_xr - pore_xl - snake_head_sigma * 4)
                                  + snake_head_sigma * 2,
                                  np.random.random() * box[1],
                                  box[2] * 0.5]

        # REENABLE EXPLOSION
        if exploding and time.time() - expl_time > 1:
            exploding = False

        # VISUALIZER
        visualizer.update()
        if has_pygame:
            if joystick_control:
                pygame.event.get()
                axis_l = np.array(
                    [joystick.get_axis(0), -joystick.get_axis(1)])
                axis_r = np.array(
                    [joystick.get_axis(3), -joystick.get_axis(4)])

                button_A = joystick.get_button(0)
                button_Start = joystick.get_button(7)

                if not button_A_old and button_A:
                    explode()
                if not button_Start_old and button_Start:
                    restart()

                button_A_old = button_A
                button_Start_old = button_A

                hat = joystick.get_hat(0)
                F_act_j = np.clip(np.array(hat) + axis_l + axis_r, -1, 1)

                set_particle_force()


t = Thread(target=main)
t.daemon = True
t.start()
visualizer.start()
