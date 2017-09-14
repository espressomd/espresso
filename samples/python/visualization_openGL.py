#!/usr/bin/python
import espressomd
from espressomd import thermostat
from espressomd import integrate
import numpy
from threading import Thread
from math import *
from espressomd.visualizationOpenGL import *

#Minimal interactive OpenGL visualization for ESPResSo

#Callbacks to control temperature 
temperature = 1.0
def increaseTemp():
        global temperature
        temperature += 0.1
        system.thermostat.set_langevin(kT=temperature, gamma=1.0)
        print "T =",system.thermostat.get_state()[0]['kT']

def decreaseTemp():
    global temperature
    temperature -= 0.1

    if temperature > 0:
        system.thermostat.set_langevin(kT=temperature, gamma=1.0)
        print "T =",system.thermostat.get_state()[0]['kT']
    else:
        temperature = 0
        system.thermostat.turn_off()
        print "T = 0"

system = espressomd.System()

system.time_step = 0.00001
system.cell_system.skin = 0.4
box_l = 10
system.box_l = [box_l, box_l, box_l]

for i in range(10):
    rpos=numpy.random.random(3) * box_l
    system.part.add(id=i, pos=rpos)
                                    
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

#system.non_bonded_inter[0,0].lennard_jones.set_params(
#    epsilon=0, sigma=1,
#    cutoff=2, shift="auto")

visualizer = openGLLive(system)

for key, value in visualizer.specs.items():
    print(str(key) + ":" + (30-len(key))*" " + str(value))

#Register buttons
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('t',KeyboardFireEvent.Hold,increaseTemp))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent('g',KeyboardFireEvent.Hold,decreaseTemp))

#Register additional callback to run in main thread
def muu():
    print "muu"
def foo():
    print "foo"

visualizer.registerCallback(foo,interval=500)
visualizer.registerCallback(muu)

def main():

    while True: 
        system.integrator.run(1)
        #Update particle information safely here
        visualizer.update()

#Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

#Start blocking visualizer
visualizer.start()
