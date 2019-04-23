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
Visualization samples with keyboard control and particle repositioning via
timed callback.
"""

from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import integrate
import numpy as np
from threading import Thread
from espressomd import visualization_opengl

required_features = []
espressomd.assert_features(required_features)

print("Press u/j to change temperature")

box_l = 10.0
system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
visualizer = visualization_opengl.openGLLive(
    system, drag_enabled=True, drag_force=100)

system.time_step = 0.00001
system.cell_system.skin = 3.0

N = 50
for i in range(N):
    system.part.add(pos=[0, 0, 0])

system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Callback for particle positions/velocities


def spin():
    system.part[:].pos = [[box_l * 0.5, box_l * (i + 1) / (N + 2), box_l * 0.5]
                          for i in range(N)]
    system.part[:].v = [
        [np.sin(10.0 * i / N) * 20, 0, np.cos(10.0 * i / N) * 20] for i in range(N)]


# Register timed callback
visualizer.register_callback(spin, interval=5000)

# Callbacks to control temperature
temperature = 1.0


def increaseTemp():
    global temperature
    temperature += 0.1
    system.thermostat.set_langevin(kT=temperature, gamma=1.0)
    print("T =", system.thermostat.get_state()[0]['kT'])


def decreaseTemp():
    global temperature
    temperature -= 0.1

    if temperature > 0:
        system.thermostat.set_langevin(kT=temperature, gamma=1.0)
        print("T =", system.thermostat.get_state()[0]['kT'])
    else:
        temperature = 0
        system.thermostat.turn_off()
        print("T = 0")


# Register button callbacks
visualizer.keyboardManager.register_button(
    visualization_opengl.KeyboardButtonEvent(
        'u', visualization_opengl.KeyboardFireEvent.Hold, increaseTemp))
visualizer.keyboardManager.register_button(
    visualization_opengl.KeyboardButtonEvent(
        'j', visualization_opengl.KeyboardFireEvent.Hold, decreaseTemp))

# Set initial position
spin()

# Start the visualizer and run the integration thread
visualizer.run(1)
