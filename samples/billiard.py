#
# Copyright (C) 2010-2022 The ESPResSo project
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

"""
ESPResSo 8Ball billiards game.
"""

import numpy as np
import math
import threading

import espressomd
import espressomd.interactions
from espressomd.visualization import openGLLive, KeyboardButtonEvent, KeyboardFireEvent
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
    impulse_step_progression = 1.2
    impulse_step_maximum = 1.5
    impulse_step_minimum = 0.5
    angle_step_progression = 1.1
    angle_step_maximum = 0.2
    angle_step_minimum = 0.02

    def __init__(self):
        self.stopped = True
        self.angle = np.pi * 0.5
        self.angle_step = 0.
        self.impulse = 10.0
        self.impulse_step = 0.

    def update_cueball_force(self):
        direction = np.array([math.sin(self.angle), 0., math.cos(self.angle)])
        cue_ball.ext_force = self.impulse * direction

    def resetAngle(self):
        self.angle_step = 0.

    def decreaseAngle(self):
        if self.stopped:
            if self.angle_step <= 0.:
                self.angle_step = self.angle_step_minimum
            elif self.angle_step < self.angle_step_maximum:
                self.angle_step *= self.angle_step_progression
            self.angle += self.angle_step
            self.update_cueball_force()

    def increaseAngle(self):
        if self.stopped:
            if self.angle_step >= 0.:
                self.angle_step = -self.angle_step_minimum
            elif self.angle_step > -self.angle_step_maximum:
                self.angle_step *= self.angle_step_progression
            self.angle += self.angle_step
            self.update_cueball_force()

    def resetImpulse(self):
        self.impulse_step = 0.

    def decreaseImpulse(self):
        if self.stopped:
            if self.impulse_step >= 0.:
                self.impulse_step = -self.impulse_step_minimum
            elif self.impulse_step > -self.impulse_step_maximum:
                self.impulse_step *= self.impulse_step_progression
            self.impulse += self.impulse_step
            self.update_cueball_force()

    def increaseImpulse(self):
        if self.stopped:
            if self.impulse_step <= 0.:
                self.impulse_step = self.impulse_step_minimum
            elif self.impulse_step < self.impulse_step_maximum:
                self.impulse_step *= self.impulse_step_progression
            self.impulse += self.impulse_step
            self.update_cueball_force()

    def fire(self):
        if self.stopped:
            self.stopped = False
            cue_ball.v = cue_ball.v + cue_ball.ext_force
            cue_ball.fix = [False, True, False]
            cue_ball.ext_force = [0, 0, 0]


pool = Billiards()

visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('4', KeyboardFireEvent.Hold, pool.decreaseAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('4', KeyboardFireEvent.Released, pool.resetAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('6', KeyboardFireEvent.Hold, pool.increaseAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('6', KeyboardFireEvent.Released, pool.resetAngle))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('2', KeyboardFireEvent.Hold, pool.decreaseImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('2', KeyboardFireEvent.Released, pool.resetImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('8', KeyboardFireEvent.Hold, pool.increaseImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('8', KeyboardFireEvent.Released, pool.resetImpulse))
visualizer.keyboard_manager.register_button(
    KeyboardButtonEvent('5', KeyboardFireEvent.Pressed, pool.fire))


system.time_step = 0.00008
system.cell_system.skin = 0.4

table_friction = 0.9
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
cue_ball = system.part.add(pos=ball_start_pos,
                           type=types['cue_ball'], mass=ball_mass)
spawnpos = {cue_ball.id: ball_start_pos}

d = 1.15 * ball_diam
a1 = np.array([d * math.sqrt(3) / 2.0, 0, -0.5 * d])
a2 = np.array([d * math.sqrt(3) / 2.0, 0, 0.5 * d])
sp = [system.box_l[0] * 0.7, ball_y,
      system.box_l[2] * 0.5 + ball_diam * 0.5]
order = (
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
    types['striped_ball'],
)

n = 0
for i in range(5):
    for j in range(i + 1):
        t = order[n]
        pos = sp + a1 * (i + 1 - j) + a2 * j
        p = system.part.add(pos=pos, mass=ball_mass,
                            type=t, fix=[False, True, False])
        spawnpos[p.id] = pos
        if t == types["black_ball"]:
            black_ball = p
        n += 1

pool.update_cueball_force()
cue_ball.fix = [True, True, True]
system.thermostat.set_langevin(kT=0., gamma=table_friction, seed=42)

cleared_balls = [0, 0]


def main():
    while True:
        system.integrator.run(1)

        vsum = 0
        for p in system.part:
            vsum += np.linalg.norm(p.v)

            for h in hole_pos:
                d = ((p.pos_folded[0] - h[0])**2
                     + (p.pos_folded[2] - h[2])**2)**0.5
                if d < hole_score_rad:
                    if p.id == cue_ball.id:
                        p.pos = ball_start_pos
                        p.v = [0, 0, 0]
                    elif p.id == black_ball.id:
                        for p2 in system.part:
                            p2.pos = spawnpos[p2.id]
                            p2.v = [0, 0, 0]
                            p2.fix = [False, True, False]
                        cue_ball.fix = [True, True, True]
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
            cue_ball.fix = [True, True, True]
            for p in system.part:
                p.v = [0, 0, 0]
            pool.update_cueball_force()

        visualizer.update()


t = threading.Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
