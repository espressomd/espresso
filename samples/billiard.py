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
ESPResSo 8Ball billiards game.
"""

import numpy as np
import math
from threading import Thread

import espressomd
import espressomd.interactions
from espressomd.visualization_opengl import openGLLive, KeyboardButtonEvent, KeyboardFireEvent
import espressomd.shapes

required_features = ["WCA", "MASS", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

print('''8Ball BILLIARDS - An ESPResSo Visualizer Demo
Controls:
  Numpad 4/6: Adjust Angle
  Numpad 2/8: Adjust Impulse
  Numpad 5: Shoot''')

# ESPRESSO
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
table_dim = [2.24, 1.12]
system.box_l = [table_dim[0], 3, table_dim[1]]

visualizer = openGLLive(
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


class Billiards:
    def __init__(self):
        self.stopped = True
        self.angle = np.pi * 0.5
        self.impulse = 10.0

    def update_cueball_force(self):
        direction = np.array([math.sin(self.angle), 0, math.cos(self.angle)])
        system.part[0].ext_force = self.impulse * direction

    def decreaseAngle(self):
        if self.stopped:
            self.angle += 0.01
            self.update_cueball_force()

    def increaseAngle(self):
        if self.stopped:
            self.angle -= 0.01
            self.update_cueball_force()

    def decreaseImpulse(self):
        if self.stopped:
            self.impulse -= 0.5
            self.update_cueball_force()

    def increaseImpulse(self):
        if self.stopped:
            self.impulse += 0.5
            self.update_cueball_force()

    def fire(self):
        if self.stopped:
            self.stopped = False
            system.part[0].v = system.part[0].v + system.part[0].ext_force
            system.part[0].fix = [False, True, False]
            system.part[0].ext_force = [0, 0, 0]


pool = Billiards()

visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('4', KeyboardFireEvent.Hold, pool.decreaseAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('6', KeyboardFireEvent.Hold, pool.increaseAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('2', KeyboardFireEvent.Hold, pool.decreaseImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('8', KeyboardFireEvent.Hold, pool.increaseImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('5', KeyboardFireEvent.Pressed, pool.fire))


def main():

    system.time_step = 0.00008
    system.cell_system.skin = 0.4

    table_h = 0.5
    ball_diam = 0.0572
    ball_mass = 0.17
    hole_dist = 0.02
    hole_rad = 0.08
    hole_score_rad = hole_rad + ball_diam / 2
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

    # WCA
    for t1 in range(4):
        for t2 in range(6):
            system.non_bonded_inter[t1, t2].wca.set_params(
                epsilon=1.0, sigma=ball_diam)

    ball_y = table_h + ball_diam * 1.5

    # PARTICLES
    ball_start_pos = [table_dim[0] * 0.25, ball_y, table_dim[1] * 0.5]
    system.part.add(id=0, pos=ball_start_pos,
                    type=types['cue_ball'], mass=ball_mass)
    spawnpos = []
    spawnpos.append(ball_start_pos)
    ball = system.part[0]

    d = 1.15 * ball_diam
    a1 = np.array([d * math.sqrt(3) / 2.0, 0, -0.5 * d])
    a2 = np.array([d * math.sqrt(3) / 2.0, 0, 0.5 * d])
    sp = [system.box_l[0] * 0.7, ball_y,
          system.box_l[2] * 0.5 + ball_diam * 0.5]
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
            system.part.add(id=pid, pos=pos, mass=ball_mass,
                            type=t, fix=[False, True, False])
            spawnpos.append(pos)
            pid += 1

    pool.update_cueball_force()
    ball.fix = [True, True, True]
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
                if d < hole_score_rad:
                    if p.id == 0:
                        p.pos = ball_start_pos
                        p.v = [0, 0, 0]
                    elif p.id == 5:
                        for p in system.part:
                            p.pos = spawnpos[p.id]
                            p.v = [0, 0, 0]
                            p.fix = [False, True, False]
                        ball.fix = [True, True, True]
                        pool.update_cueball_force()
                        pool.stopped = True
                    else:
                        t = p.type - 1
                        cleared_balls[t] += 1
                        if t == 0:
                            z = table_dim[1] - ball_diam * 0.6
                        else:
                            z = ball_diam * 0.6
                        p.pos = [cleared_balls[t] * ball_diam * 1.5, 1.1, z]
                        p.fix = [True, True, True]
                        p.v = [0, 0, 0]

        if not pool.stopped and vsum < 0.3:
            pool.stopped = True
            ball.fix = [True, True, True]
            for p in system.part:
                p.v = [0, 0, 0]
            pool.update_cueball_force()

        visualizer.update()


t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
