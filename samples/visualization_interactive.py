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
Visualize a simulation box where the particles can be repositioned via the
mouse and timed callbacks, and the temperature of the thermostat can be
changed via the keyboard.
"""

import espressomd
import espressomd.visualization_opengl
import numpy as np

required_features = []
espressomd.assert_features(required_features)

print("Press u/j to change temperature")

box_l = 10.0
system = espressomd.System(box_l=[box_l] * 3)
visualizer = espressomd.visualization_opengl.openGLLive(
    system, drag_enabled=True, drag_force=100)

system.time_step = 0.00001
system.cell_system.skin = 3.0

N = 50
partcls = system.part.add(pos=N * [[0, 0, 0]])

system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Callback for particle positions/velocities


def spin():
    partcls.pos = [[box_l * 0.5, box_l * (i + 1) / (N + 2), box_l * 0.5]
                   for i in range(N)]
    partcls.v = [
        [np.sin(10.0 * i / N) * 20, 0, np.cos(10.0 * i / N) * 20] for i in range(N)]


# Register timed callback
visualizer.register_callback(spin, interval=5000)

# Callbacks to control temperature
temperature = 1.0


def increaseTemp():
    global temperature
    temperature += 0.5
    system.thermostat.set_langevin(kT=temperature, gamma=1.0)
    print("T =", system.thermostat.get_state()[0]['kT'])


def decreaseTemp():
    global temperature
    temperature -= 0.5

    if temperature > 0:
        system.thermostat.set_langevin(kT=temperature, gamma=1.0)
        print(f"T = {system.thermostat.get_state()[0]['kT']:.1f}")
    else:
        temperature = 0
        system.thermostat.turn_off()
        print("T = 0.")


# Register button callbacks
visualizer.keyboard_manager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent('u', espressomd.visualization_opengl.KeyboardFireEvent.Hold, increaseTemp))
visualizer.keyboard_manager.register_button(
    espressomd.visualization_opengl.KeyboardButtonEvent('j', espressomd.visualization_opengl.KeyboardFireEvent.Hold, decreaseTemp))

# Set initial position
spin()

# Start the visualizer
visualizer.run(1)
