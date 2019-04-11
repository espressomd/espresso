# Copyright (C) 2010-2018 The ESPResSo project
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
ESPResSo 8Ball billiard game.
"""

from __future__ import print_function
import numpy as np
import math
from threading import Thread

import espressomd
from espressomd import thermostat
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
from espressomd import minimize_energy
import espressomd.interactions
import espressomd.visualization_opengl
import espressomd.shapes

required_features = ["LENNARD_JONES", "MASS", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

print('''8Ball BILLIARD - An ESPResSo Visualizer Demo
Controls:
  Numpad 4/6: Adjust Angle
  Numpad 2/8: Adjust Impulse
  Numpad 5: Shoot''')

# ESPRESSO
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
table_dim = [2.24, 1.12]
system.box_l = [table_dim[0], 3, table_dim[1]]

visualizer = espressomd.visualization_opengl.openGLLive(
    system,
    ext_force_arrows=True,
    ext_force_arrows_type_scale=[0.02],
    ext_force_arrows_type_radii=[0.01],
    background_color=[0.5, 0.4, 0.5],
    drag_enabled=False,
    particle_type_materials=['medium', 'bright', 'bright', 'medium'],
    particle_type_colors=[
        [1, 1, 1], [0.5, 0.1, 0.1], [0.1, 0.2, 0.4], [0.2, 0.2, 0.2]],
    constraint_type_materials=['dark'],
    constraint_type_colors=[[0.1, 0.424, 0.011], [0.1, 0.1, 0.1]],
    camera_position=[1.12, 2.8, 0.56],
    window_size=[1000, 600],
    draw_axis=False,
    light_pos=[table_dim[0] * 0.5, 1.0, table_dim[1] * 0.5],
    light_colors=[[0.8, 0.8, 0.8], [0.9, 0.9, 0.9], [1.0, 1.0, 1.0]],
    light_brightness=1.0)

stopped = True
angle = np.pi * 0.5
impulse = 10.0


def decreaseAngle():
    global angle, impulse
    if stopped:
        angle += 0.01
        system.part[0].ext_force = impulse * \
            np.array([math.sin(angle), 0, math.cos(angle)])


def increaseAngle():
    global angle, impulse
    if stopped:
        angle -= 0.01
        system.part[0].ext_force = impulse * \
            np.array([math.sin(angle), 0, math.cos(angle)])


def decreaseImpulse():
    global impulse, angle
    if stopped:
        impulse -= 0.5
        system.part[0].ext_force = impulse * \
            np.array([math.sin(angle), 0, math.cos(angle)])


def increaseImpulse():
    global impulse, angle
    if stopped:
        impulse += 0.5
        system.part[0].ext_force = impulse * \
            np.array([math.sin(angle), 0, math.cos(angle)])


def fire():
    global stopped
    if stopped:
        stopped = False
        system.part[0].v = system.part[0].v + \
            impulse * np.array([math.sin(angle), 0, math.cos(angle)])
        system.part[0].fix = [0, 1, 0]
        system.part[0].ext_force = [0, 0, 0]


visualizer.keyboardManager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent(
        '4', espressomd.visualization_opengl.KeyboardFireEvent.Hold, decreaseAngle))
visualizer.keyboardManager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent(
        '6', espressomd.visualization_opengl.KeyboardFireEvent.Hold, increaseAngle))
visualizer.keyboardManager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent(
        '2', espressomd.visualization_opengl.KeyboardFireEvent.Hold, decreaseImpulse))
visualizer.keyboardManager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent(
        '8', espressomd.visualization_opengl.KeyboardFireEvent.Hold, increaseImpulse))
visualizer.keyboardManager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent(
        '5', espressomd.visualization_opengl.KeyboardFireEvent.Pressed, fire))


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
                [table_dim[0] * 0.5, table_h, table_dim[1] - hole_dist],
                [table_dim[0] * 0.5, table_h, hole_dist]]
    types = {'cue_ball': 0, 'striped_ball': 1, 'solid_ball': 2,
             'black_ball': 3, 'table': 4, 'wall': 5, 'hole': 6}

    system.constraints.add(
        shape=espressomd.shapes.Wall(dist=table_h, normal=[0.0, 1.0, 0.0]),
        particle_type=types['table'],
        penetrable=True)

    system.constraints.add(
        shape=espressomd.shapes.Wall(dist=0.01, normal=[1.0, 0.0, 0.0]),
        particle_type=types['wall'],
        penetrable=True)
    system.constraints.add(
        shape=espressomd.shapes.Wall(
            dist=-(table_dim[0] - 0.01),
            normal=[-1.0, 0.0, 0.0]),
        particle_type=types['wall'],
        penetrable=True)
    system.constraints.add(
        shape=espressomd.shapes.Wall(
            dist=0.01,
            normal=[0.0, 0.0, 1.0]),
        particle_type=types['wall'],
        penetrable=True)
    system.constraints.add(
        shape=espressomd.shapes.Wall(
            dist=-(table_dim[1] - 0.01),
            normal=[0.0, 0.0, -1.0]),
        particle_type=types['wall'],
        penetrable=True)
    for h in hole_pos:
        system.constraints.add(
            shape=espressomd.shapes.Cylinder(
                center=(np.array(h)
                        - np.array([0, table_h * 0.5, 0])).tolist(),
                axis=[0, 1, 0],
                radius=hole_rad,
                length=1.02 * table_h,
                direction=1),
            particle_type=types['hole'],
            penetrable=True)

    lj_eps = np.array([1])
    lj_sig = np.array([ball_diam])
    lj_cut = lj_sig * 2.0**(1.0 / 6.0)
    lj_cap = 20
    mass = np.array([0.17])

    num_types = len(lj_sig)

    # LENNARD JONES
    def mix_eps(eps1, eps2, rule='LB'):
        return math.sqrt(eps1 * eps2)

    def mix_sig(sig1, sig2, rule='LB'):
        return 0.5 * (sig1 + sig2)

    for t1 in range(4):
        for t2 in range(6):
            system.non_bonded_inter[t1, t2].lennard_jones.set_params(
                epsilon=mix_eps(lj_eps[0], lj_eps[0]),
                sigma=mix_sig(lj_sig[0], lj_sig[0]),
                cutoff=mix_sig(lj_cut[0], lj_cut[0]),
                shift="auto")

    ball_y = table_h + ball_diam * 1.5

    # PARTICLES
    ball_start_pos = [table_dim[0] * 0.25, ball_y, table_dim[1] * 0.5]
    system.part.add(id=0, pos=ball_start_pos,
                    type=types['cue_ball'], mass=mass[0])
    spawnpos = []
    spawnpos.append(ball_start_pos)
    ball = system.part[0]

    d = lj_sig[0] * 1.15
    a1 = np.array([d * math.sqrt(3) / 2.0, 0, -0.5 * d])
    a2 = np.array([d * math.sqrt(3) / 2.0, 0, 0.5 * d])
    sp = [system.box_l[0] * 0.7, ball_y,
          system.box_l[2] * 0.5 + lj_sig[0] * 0.5]
    pid = 1
    order = [
        types['solid_ball'],
        types['striped_ball'],
        types['solid_ball'],
        types['solid_ball'],
        types['black_ball'],
        types['striped_ball'],
        types['striped_ball'],
        types['solid_ball'],
        types['striped_ball'],
        types['solid_ball'],
        types['solid_ball'],
        types['striped_ball'],
        types['striped_ball'],
        types['solid_ball'],
        types['striped_ball']]

    for i in range(5):
        for j in range(i + 1):
            N = i + 1
            t = order[pid - 1]
            pos = sp + a1 * (N - j) + a2 * j
            system.part.add(
                id=pid, pos=pos, mass=mass[0], type=t, fix=[0, 1, 0])
            spawnpos.append(pos)
            pid += 1

    ball.ext_force = impulse * np.array([math.sin(angle), 0, math.cos(angle)])
    ball.fix = [1, 1, 1]
    system.thermostat.set_langevin(kT=0, gamma=0.8, seed=42)

    cleared_balls = [0, 0]
    while True:
        system.integrator.run(1)

        vsum = 0
        for p in system.part:
            vsum += np.linalg.norm(p.v)

            for h in hole_pos:

                d = ((p.pos_folded[0] - h[0])**2
                     + (p.pos_folded[2] - h[2])**2)**0.5
                if (d < hole_score_rad):
                    if p.id == 0:
                        p.pos = ball_start_pos
                        p.v = [0, 0, 0]
                    elif p.id == 5:
                        for p in system.part:
                            p.pos = spawnpos[p.id]
                            p.v = [0, 0, 0]
                            p.fix = [0, 1, 0]
                        ball.fix = [1, 1, 1]
                        ball.ext_force = impulse * \
                            np.array([math.sin(angle), 0, math.cos(angle)])
                        stoppen = True
                    else:
                        t = p.type - 1
                        cleared_balls[t] += 1
                        if t == 0:
                            z = table_dim[1] - lj_sig[0] * 0.6
                        else:
                            z = lj_sig[0] * 0.6
                        p.pos = [cleared_balls[t] * lj_sig[0] * 1.5, 1.1, z]
                        p.fix = [1, 1, 1]
                        p.v = [0, 0, 0]

        if not stopped and vsum < 0.3:
            stopped = True
            ball.fix = [1, 1, 1]
            for p in system.part:
                p.v = [0, 0, 0]
            ball.ext_force = impulse * \
                np.array([math.sin(angle), 0, math.cos(angle)])

        visualizer.update()


t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
