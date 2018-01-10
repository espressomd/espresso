from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import integrate
import numpy
from threading import Thread
from math import *
from espressomd.visualization_opengl import *

# Minimal interactive OpenGL visualization for ESPResSo

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


system = espressomd.System()

system.time_step = 0.00001
system.cell_system.skin = 0.4
box_l = 10
system.box_l = [box_l, box_l, box_l]

for i in range(10):
    rpos = numpy.random.random(3) * box_l
    system.part.add(id=i, pos=rpos)

system.thermostat.set_langevin(kT=1.0, gamma=1.0)

visualizer = openGLLive(system, drag_enabled = True)

for key, value in visualizer.specs.items():
    print(str(key) + ":" + (30 - len(key)) * " " + str(value))

# Register buttons
visualizer.keyboardManager.register_button(
    KeyboardButtonEvent('u', KeyboardFireEvent.Hold, increaseTemp))
visualizer.keyboardManager.register_button(
    KeyboardButtonEvent('j', KeyboardFireEvent.Hold, decreaseTemp))


# Register additional callback to run in main thread
def muu():
    print("muu")


def foo():
    print("foo")


visualizer.register_callback(foo, interval=500)
visualizer.register_callback(muu)


def main():

    while True:
        system.integrator.run(1)
        # Update particle information safely here
        visualizer.update()


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

# Start blocking visualizer
visualizer.start()
