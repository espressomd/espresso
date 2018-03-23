from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import integrate
import numpy as np
from threading import Thread
from math import *
from espressomd.visualization_opengl import *

# Minimal interactive OpenGL visualization for ESPResSo

print("Press u/j to change temperature")

box_l = 10.0
system = espressomd.System(box_l=[box_l] * 3)
visualizer = openGLLive(system, drag_enabled=True, drag_force=100)

system.time_step = 0.00001
system.cell_system.skin = 3.0

N = 50
for i in range(N):
    system.part.add(pos=[0, 0, 0])

system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# Callback for particle positions/velocities


def spin():
    system.part[:].pos = [[box_l * 0.5, box_l *
                           (i + 1) / (N + 2), box_l * 0.5] for i in range(N)]
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
    KeyboardButtonEvent('u', KeyboardFireEvent.Hold, increaseTemp))
visualizer.keyboardManager.register_button(
    KeyboardButtonEvent('j', KeyboardFireEvent.Hold, decreaseTemp))

# Set initial position
spin()

# Start the visualizer and run the integration thread
visualizer.run(1)
