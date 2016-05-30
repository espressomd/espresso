#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from __future__ import print_function
import espressomd._system as es
import espressomd
from espressomd import thermostat
from espressomd import analyze
from espressomd import integrate
from espressomd import visualization
import numpy
from matplotlib import pyplot
from threading import Thread
from traits.api import HasTraits, Button, Any, Range, List
from traitsui.api import View, Group, Item, CheckListEditor, RangeEditor
from pygame import midi

midi.init()

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 10
int_n_times = 5000000


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

analyze.distto(system, 0)

act_min_dist = analyze.mindist(system)
system.max_num_cells = 2744

mayavi = visualization.mayavi_live(system)

#############################################################
#  GUI Controls                                            #
#############################################################

inputs, outputs = [],[]
for i in range(midi.get_count()):
	interf, name, input, output, opened = midi.get_device_info(i)
	if input:
		inputs.append((i, interf + " " + name))
	if output:
		outputs.append((i, interf + " " + name))

class Controls(HasTraits):
	input_device = List(editor=CheckListEditor(values=inputs))
	output_device = List(editor=CheckListEditor(values=outputs))
	
	max_temp  = 2.
	max_press = 2.
	max_vol   = 5.
	min_vol   = 1.
	max_n     = 5000
	
	temperature = Range(0., max_temp, 1., )
	volume = Range(min_vol, max_vol, 1., )
	pressure = Range(0., max_press, 1., )
	number_of_particles = Range(1, max_n, n_part, )
	
	midi_input = midi.Input(inputs[0][0]) if len(inputs) > 0 else None
	midi_output = midi.Output(outputs[0][0]) if len(outputs) > 0 else None
	
	MIDI_CC   = 176
	MIDI_BASE = 81
	MIDI_NUM_TEMPERATURE       = MIDI_BASE+0
	MIDI_NUM_VOLUME            = MIDI_BASE+1
	MIDI_NUM_PRESSURE          = MIDI_BASE+2
	MIDI_NUM_NUMBEROFPARTICLES = MIDI_BASE+3
	
	_ui = Any
	view = View(
		Group(
			Item('temperature', editor=RangeEditor(high_name='max_temp')),
			Item('volume', editor=RangeEditor(low_name='min_vol', high_name='max_vol')),
#			Item('pressure', editor=RangeEditor(high_name='max_press')),
			Item('number_of_particles', editor=RangeEditor(high_name='max_n', is_float=False)),
			show_labels=True,
			label='Parameters'
		),
		Group(
			Item('input_device'),
			Item('output_device'),
			show_labels=True,
			label='MIDI devices'
		),
		buttons=[],
		title='Control'
	)
	
	def __init__(self, **traits):
		super(Controls, self).__init__(**traits)
		self._ui = self.edit_traits()
		
		# send the default values to the MIDI controller
		self._temperature_fired()
		self._volume_fired()
		self._pressure_fired()
		self._number_of_particles_fired()
	
	def _input_device_fired(self):
		if self.midi_input is not None:
			self.midi_input.close()
		self.midi_input = midi.Input(self.input_device[0])
	
	def _output_device_fired(self):
		if self.midi_input is not None:
			self.midi_output.close()
		self.midi_output = midi.Output(self.output_device[0])
	
	def _temperature_fired(self):
		status = self.MIDI_CC
		data1  = self.MIDI_NUM_TEMPERATURE
		data2  = int(self.temperature/self.max_temp*127)
		if self.midi_output is not None:
			self.midi_output.write_short(status,data1,data2)
	
	def _volume_fired(self):
		status = self.MIDI_CC
		data1  = self.MIDI_NUM_VOLUME
		data2  = int(self.volume/self.max_vol*127)
		if self.midi_output is not None:
			self.midi_output.write_short(status,data1,data2)
	
	def _pressure_fired(self):
		status = self.MIDI_CC
		data1  = self.MIDI_NUM_PRESSURE
		data2  = int(self.pressure/self.max_press*127)
		if self.midi_output is not None:
			self.midi_output.write_short(status,data1,data2)
	
	def _number_of_particles_fired(self):
		status = self.MIDI_CC
		data1  = self.MIDI_NUM_NUMBEROFPARTICLES
		data2  = int(self.number_of_particles/self.max_n*127)
		if self.midi_output is not None:
			self.midi_output.write_short(status,data1,data2)

#############################################################
#  Warmup Integration                                       #
#############################################################

# set LJ cap
lj_cap = 20
system.non_bonded_inter.set_force_cap(lj_cap)

# # Warmup Integration Loop
# i = 0
# while (i < warm_n_times and act_min_dist < min_dist):
#     integrate.integrate(warm_steps)
#     # Warmup criterion
#     act_min_dist = analyze.mindist(system)
#     i += 1
#
# #   Increase LJ cap
#     lj_cap = lj_cap + 10
#     system.non_bonded_inter.set_force_cap(lj_cap)
#     mayavi.update()

#############################################################
#      Integration                                          #
#############################################################

# remove force capping
#lj_cap = 0
#system.non_bonded_inter.set_force_cap(lj_cap)

# get initial observables
pressure = analyze.pressure(system)
temperature = system.temperature

# TODO: this is some terrible polynomial fit, replace it with a better expression
# equation of state
pyplot.subplot(131)
pyplot.semilogy()
pyplot.title("Phase diagram")
pyplot.xlabel("Temperature")
pyplot.ylabel("Pressure")
pyplot.xlim(0.5,2.0)
pyplot.ylim(5e-5,2e1)
xx = numpy.linspace(0.5,0.7,200)
pyplot.plot(xx, -6.726* xx**4 + 16.92* xx**3 - 15.85* xx**2 + 6.563* xx - 1.015, 'k-')
xx = numpy.linspace(0.7,1.3,600)
pyplot.plot(xx, -0.5002* xx**4 + 2.233* xx**3 - 3.207* xx**2 + 1.917* xx - 0.4151, 'k-')
xx = numpy.linspace(0.6,2.2,1500)
pyplot.plot(xx, 16.72* xx**4 - 88.28* xx**3 + 168* xx**2 - 122.4* xx +29.79, 'k-')

cursor = pyplot.scatter(temperature,pressure['total'],200,'r')
pyplot.text(0.6, 10, 'solid')
pyplot.text(1, 1, 'liquid')
pyplot.text(1, 10**-3, 'gas')

pyplot.subplot(132)
pyplot.title("Temperature")
plot1, = pyplot.plot([0],[temperature])
pyplot.xlabel("Time")
pyplot.ylabel("Temperature")
pyplot.subplot(133)
pyplot.title("Pressure")
plot2, = pyplot.plot([0],[pressure['total']])
pyplot.xlabel("Time")
pyplot.ylabel("Pressure")
#pyplot.legend()
pyplot.show(block=False)

def main_loop():
	global energies

	integrate.integrate(int_steps)
	mayavi.update()
	
	# update the parameters set in the GUI
	system.thermostat.set_langevin(kT=controls.temperature, gamma=1.0)
	
	# make sure the parameters are valid
	if controls.volume == 0:
		controls.volume = 1
	if controls.number_of_particles == 0:
		controls.number_of_particles = 1
	
	new_box = numpy.ones(3) * box_l * controls.volume
	if numpy.any(numpy.array(system.box_l) != new_box):
		for i in range(system.n_part):
			system.part[i].pos *= new_box/system.box_l
	system.box_l = new_box
	
	new_part = controls.number_of_particles
	if new_part > system.n_part:
		for i in range(system.n_part, new_part):
			system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)
	elif new_part < system.n_part:
		for i in range(new_part, system.n_part):
			system.part[i].delete()
	assert system.n_part == system.max_part +1 # There should be no gaps in particle numbers

	energies = analyze.energy(system=system)
	plot1.set_xdata(numpy.append(plot1.get_xdata(), system.time))
	plot1.set_ydata(numpy.append(plot1.get_ydata(), system.temperature))
	plot2.set_xdata(numpy.append(plot2.get_xdata(), system.time))
	plot2.set_ydata(numpy.append(plot2.get_ydata(), analyze.pressure(system)['total']))
	cursor.set_offsets(([[plot1.get_ydata()[-1], plot2.get_ydata()[-1]]]))

def main_thread():
    for i in range(0, int_n_times):
        main_loop()

def midi_thread():
	while True:
		try:
			if controls.midi_input is not None and controls.midi_input.poll():
				events = controls.midi_input.read(1000)
				for event in events:
					status,data1,data2,data3 = event[0]
					if status == controls.MIDI_CC:
						if data1 == controls.MIDI_NUM_TEMPERATURE:
							temperature = data2 * controls.max_temp / 127
							controls.temperature = temperature
						elif data1 == controls.MIDI_NUM_VOLUME:
							volume = data2 * controls.max_vol / 127
							controls.volume = volume
						elif data1 == controls.MIDI_NUM_PRESSURE:
							pressure = data2 * controls.max_press / 127
							controls.pressure = pressure
						elif data1 == controls.MIDI_NUM_NUMBEROFPARTICLES:
							npart = data2 * controls.max_n / 127
							controls.number_of_particles = npart
		except Exception as e:
			print (e)

last_plotted = 0
def update_plot():
    global last_plotted
    current_time = plot1.get_xdata()[-1]
    if last_plotted == current_time:
        return
    last_plotted = current_time
    plot1.axes.set_xlim(0, plot1.get_xdata()[-1])
    plot1.axes.set_ylim(0.8*plot1.get_ydata().min(), 1.2*plot1.get_ydata().max())
    plot2.axes.set_xlim(0, plot2.get_xdata()[-1])
    plot2.axes.set_ylim(0.8*plot2.get_ydata().min(), 1.2*plot2.get_ydata().max())
    pyplot.draw()

t = Thread(target=main_thread)
t.daemon = True
mayavi.register_callback(update_plot, interval=2000)
controls = Controls()
t.start()
if controls.midi_input is not None:
	t2 = Thread(target=midi_thread)
	t2.daemon = True
	t2.start()
mayavi.run_gui_event_loop()
