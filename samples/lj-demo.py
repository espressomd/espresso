#
# Copyright (C) 2013-2019 The ESPResSo project
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
Simulate a Lennard-Jones fluid in different thermodynamic ensembles (NVT, NpT).
Sliders from a MIDI controller can change system variables such as temperature
and volume. Some thermodynamic observables are analyzed and plotted live.
"""
import matplotlib
matplotlib.use('WXAgg')
import espressomd
espressomd.assert_features(["LENNARD_JONES"])
from espressomd import visualization
import numpy as np
from matplotlib import pyplot
from threading import Thread
from traits.api import HasTraits, Any, Range, List, Enum, Float
from traitsui.api import View, Group, Item, CheckListEditor, RangeEditor
import time
import argparse

parser = argparse.ArgumentParser(epilog=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument("--mayavi", action="store_const", dest="visualizer",
                   const="mayavi", help="MayaVi visualizer", default="mayavi")
group.add_argument("--opengl", action="store_const", dest="visualizer",
                   const="opengl", help="OpenGL visualizer")
args = parser.parse_args()

use_opengl = args.visualizer == "opengl"
use_mayavi = args.visualizer == "mayavi"

if use_mayavi:
    from espressomd.visualization_mayavi import mlab
if use_opengl:
    from pyface.api import GUI

try:
    import midi
except BaseException:
    try:
        from pygame import midi
    except BaseException:
        from portmidi import midi
midi.init()

# if log flag is set, midi controller will change pressure logarithmically
pressure_log_flag = True

mayavi_autozoom = False  # autozoom is buggy... works only for rotation

old_pressure = -1


# NPT variables
#############################################################
NPTGamma0 = 1.0
#NPTInitPistonMass = 1e-06
#NPTMinPistonMass =  1e-06
NPTMinPistonMass = 1e-04
NPTMaxPistonMass = 1.0
NPTInitPistonMass = NPTMinPistonMass


# System parameters
#############################################################

# 300  Particles
box_l = 7.5395
density = 0.7

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5 * lj_sig
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l, box_l, box_l])

system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

system.cell_system.set_n_square(use_verlet_lists=False)

# do the warmup until the particles have at least the distance min_dist
min_dist = 0.9

# integration
int_steps = 1
int_n_times = 5000000


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

# Particle setup
#############################################################

volume = box_l**3
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

act_min_dist = system.analysis.min_dist()

if use_mayavi:
    vis = visualization.mayaviLive(system)
elif use_opengl:
    vis = visualization.openGLLive(system)

mayavi_rotation_angle = 45.
mayavi_rotation_angle_step = 5.

mayavi_zoom = 36.
mayavi_zoom_old = mayavi_zoom
mayavi_zoom_step = 3.

plot_max_data_len = 20

#############################################################
#  GUI Controls                                            #
#############################################################

inputs, outputs = [], []
for i in range(midi.get_count()):
    interf, name, input, output, opened = midi.get_device_info(i)
    if input:
        inputs.append((i, interf + " " + name))
    if output:
        outputs.append((i, interf + " " + name))


class Controls(HasTraits):
    if len(inputs) == 1:
        default_input = inputs

    for i in inputs:
        if "Through Port" not in i[1]:
            default_input = i
            break

    default_input = default_input if inputs else None

    default_output = -1
    through_port_output = None
    for i in outputs:
        if "Through Port" not in i[1]:
            default_output = i
            break
        else:
            through_port_output = i
    default_output = default_output if len(
        outputs) > 1 else through_port_output

    if default_input is None or default_output is None:
        print('Cannot connect to any MIDI device')

    input_device = List(value=default_input,
                        editor=CheckListEditor(values=inputs))
    output_device = List(value=default_output,
                         editor=CheckListEditor(values=outputs))

    max_temp = 2.
    min_temp = 0.5
    max_press = 10.
    min_press = 5e-4
    max_vol = 100000.
    min_vol = 50.
    max_n = 1000
    min_n = 50

    temperature = Range(min_temp, max_temp, 1., )
    volume = Float(box_l**3.)
    pressure = Float(1.)
    number_of_particles = Range(min_n, max_n, n_part, )
    ensemble = Enum('NVT', 'NPT')

    midi_input = None
    midi_output = None

    MIDI_BASE = 224
    MIDI_NUM_TEMPERATURE = MIDI_BASE + 0
    MIDI_NUM_VOLUME = MIDI_BASE + 1
    MIDI_NUM_PRESSURE = MIDI_BASE + 2
    MIDI_NUM_NUMBEROFPARTICLES = MIDI_BASE + 3

    MIDI_ROTATE = 0

    MIDI_ZOOM = 144

    _ui = Any
    view = View(
        Group(
            Item('temperature', editor=RangeEditor(
                low_name='min_temp', high_name='max_temp')),
            Item('volume', editor=RangeEditor(
                low_name='min_vol', high_name='max_vol')),
            Item('pressure', editor=RangeEditor(
                low_name='min_press', high_name='max_press')),
            Item('number_of_particles', editor=RangeEditor(
                low_name='min_n', high_name='max_n', is_float=False)),
            Item('ensemble', style='custom'),
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
        title='Control',
        height=0.2,
        width=0.3
    )

    def __init__(self, **traits):
        super(Controls, self).__init__(**traits)
        self._ui = self.edit_traits()
        self.push_current_values()

    def push_current_values(self):
        """send the current values to the MIDI controller"""
        self._temperature_fired()
        self._volume_fired()
        self._pressure_fired()
        self._number_of_particles_fired()
        self._ensemble_fired()

    def _input_device_fired(self):
        if self.midi_input is not None:
            self.midi_input.close()
        if self.input_device:
            self.midi_input = midi.Input(self.input_device[0])

    def _output_device_fired(self):
        if self.midi_output is not None:
            self.midi_output.close()
        self.midi_output = midi.Output(self.output_device[0])
        self.push_current_values()

    def _temperature_fired(self):
        status = self.MIDI_NUM_TEMPERATURE
        data1 = int((self.temperature - self.min_temp) /
                    (self.max_temp - self.min_temp) * 127)
        data2 = data1
        if self.midi_output is not None:
            self.midi_output.write_short(status, data1, data2)

    def _volume_fired(self):
        status = self.MIDI_NUM_VOLUME
        data1 = limit_range(int((system.box_l[0]**3. - self.min_vol) / (
            self.max_vol - self.min_vol) * 127), minval=0, maxval=127)
        data2 = data1

        if self.midi_output is not None:
            self.midi_output.write_short(status, data1, data2)

    def _pressure_fired(self):
        status = self.MIDI_NUM_PRESSURE

        if pressure_log_flag:
            data1 = limit_range(int(127 *
                                    (np.log(self.pressure) -
                                     np.log(self.min_press)) /
                                    (np.log(self.max_press) -
                                        np.log(self.min_press))), minval=0, maxval=127)
        else:
            data1 = limit_range(int((self.pressure -
                                     self.min_press) /
                                    (self.max_press -
                                     self.min_press) *
                                    127), minval=0, maxval=127)
        data2 = data1
        if self.midi_output is not None:
            self.midi_output.write_short(status, data1, data2)

    def _number_of_particles_fired(self):
        status = self.MIDI_NUM_NUMBEROFPARTICLES
        data1 = int(self.number_of_particles / self.max_n * 127)
        data2 = data1
        if self.midi_output is not None:
            self.midi_output.write_short(status, data1, data2)

    def _ensemble_fired(self):
        if self.midi_output is not None:
            self.midi_output.write_short(144, 0, 127)  # T
            self.midi_output.write_short(
                144, 1, 127 * (self.ensemble != 'NPT'))  # V
            self.midi_output.write_short(
                144, 2, 127 * (self.ensemble == 'NPT'))  # P
            self.midi_output.write_short(144, 3, 127)  # N

#############################################################
#      Integration                                          #
#############################################################


# get initial observables
pressure = system.analysis.pressure()
temperature = 0.0

# TODO: this is some terrible polynomial fit, replace it with a better expression
# equation of state
pyplot.subplot(131)
pyplot.semilogy()
pyplot.title("Phase diagram")
pyplot.xlabel("Temperature")
pyplot.ylabel("Pressure")
pyplot.xlim(0.5, 2.0)
pyplot.ylim(5e-5, 2e1)
xx = np.linspace(0.5, 0.7, 200)
pyplot.plot(xx, -6.726 * xx**4 + 16.92 * xx**3 -
            15.85 * xx**2 + 6.563 * xx - 1.015, 'k-')
xx = np.linspace(0.7, 1.3, 600)
pyplot.plot(xx, -0.5002 * xx**4 + 2.233 * xx**3 -
            3.207 * xx**2 + 1.917 * xx - 0.4151, 'k-')
xx = np.linspace(0.6, 2.2, 1500)
pyplot.plot(xx, 16.72 * xx**4 - 88.28 * xx**3 +
            168 * xx**2 - 122.4 * xx + 29.79, 'k-')

cursor = pyplot.scatter(temperature, pressure['total'], 200, 'g')
#cursor2 = pyplot.scatter(-1, -1, 200, 'r')
pyplot.text(0.6, 10, 'solid')
pyplot.text(1, 1, 'liquid')
pyplot.text(1, 10**-3, 'gas')

pyplot.subplot(132)
pyplot.title("Temperature")
plot1, = pyplot.plot([0], [temperature])
pyplot.xlabel("Time")
pyplot.ylabel("Temperature")
pyplot.subplot(133)
pyplot.title("Pressure")
plot2, = pyplot.plot([0], [pressure['total']])
pyplot.xlabel("Time")
pyplot.ylabel("Pressure")
# pyplot.legend()
pyplot.show(block=False)

plt1_x_data = np.zeros(1)
plt1_y_data = np.zeros(1)
plt2_x_data = np.zeros(1)
plt2_y_data = np.zeros(1)


def limit_range(val, minval=0., maxval=1.):
    if val > maxval:
        ret_val = maxval
    elif val < minval:
        ret_val = minval
    else:
        ret_val = val

    if isinstance(val, int):
        return int(ret_val)
    elif isinstance(val, float):
        return float(ret_val)
    else:
        return ret_val


def pressure_from_midi_val(midi_val, pmin, pmax, log_flag=pressure_log_flag):
    if log_flag:
        return pmin * (float(pmax) / pmin)**(float(midi_val) / 127)
    else:
        return midi_val * (pmax - pmin) / 127 + pmin


def main_loop():
    global energies, plt1_x_data, plt1_y_data, plt2_x_data, plt2_y_data, old_pressure

    system.integrator.run(steps=int_steps)
    vis.update()

    # increase LJ cap during warmup
    if system.force_cap > 0:
        if system.analysis.min_dist() < min_dist:
            system.force_cap = system.force_cap + 0.1
        else:
            system.force_cap = 0
            print("Switching off force capping")

    # make sure the parameters are valid
    # not sure if this is necessary after using limit_range
    if controls.volume == 0:
        controls.volume = controls.min_vol
    if controls.number_of_particles == 0:
        controls.number_of_particles = 1
    if controls.pressure == 0:
        controls.pressure = controls.min_press

    pressure = system.analysis.pressure()

    # update the parameters set in the GUI
    if system.thermostat.get_state()[0]['kT'] != controls.temperature:
        system.thermostat.set_langevin(kT=controls.temperature, gamma=1.0)
        print("temperature changed")
        system.force_cap = lj_cap
    if controls.ensemble == 'NPT':
        # reset Vkappa when target pressure has changed

        if old_pressure != controls.pressure:
            system.analysis.v_kappa('reset')
            print("pressure changed")
            old_pressure = controls.pressure
            system.force_cap = lj_cap

        newVkappa = system.analysis.v_kappa('read')['Vk1']
        newVkappa = newVkappa if newVkappa > 0. else 4.0 / \
            (NPTGamma0 * NPTGamma0 * NPTInitPistonMass)
        pistonMass = limit_range(4.0 / (NPTGamma0 * NPTGamma0 * newVkappa),
                                 NPTMinPistonMass, NPTMaxPistonMass)
        system.integrator.set_isotropic_npt(
            controls.pressure, pistonMass, cubic_box=True)

        controls.volume = system.box_l[0]**3.

    else:
        system.integrator.set_nvt()
        controls.pressure = pressure['total']

        new_box = np.ones(3) * controls.volume**(1. / 3.)
        if np.any(np.array(system.box_l) != new_box):
            for i in range(len(system.part)):
                system.part[i].pos = system.part[i].pos * \
                    new_box / system.box_l[0]
            print("volume changed")
            system.force_cap = lj_cap
        system.box_l = new_box

    new_part = controls.number_of_particles
    if new_part > len(system.part):
        for i in range(len(system.part), new_part):
            system.part.add(id=i, pos=np.random.random(3) * system.box_l)
        print("particles added")
        system.force_cap = lj_cap
    elif new_part < len(system.part):
        for i in range(new_part, len(system.part)):
            system.part[i].remove()
        print("particles removed")

    plt1_x_data = plot1.get_xdata()
    plt1_y_data = plot1.get_ydata()
    plt2_x_data = plot2.get_xdata()
    plt2_y_data = plot2.get_ydata()

    plt1_x_data = np.append(
        plt1_x_data[-plot_max_data_len + 1:], system.time)
    plt1_y_data = np.append(plt1_y_data[-plot_max_data_len + 1:],
                            2. / (3. * len(system.part))
                            * system.analysis.energy()["kinetic"])
    plt2_x_data = np.append(
        plt2_x_data[-plot_max_data_len + 1:], system.time)
    plt2_y_data = np.append(
        plt2_y_data[-plot_max_data_len + 1:], pressure['total'])


def main_thread():
    for _ in range(int_n_times):
        main_loop()


def midi_thread():
    global mayavi_rotation_angle, mayavi_zoom
    while True:
        try:
            if controls.midi_input is not None and controls.midi_input.poll():
                events = controls.midi_input.read(1000)
                for event in events:
                    status, data1, data2, _ = event[0]
                    if status == controls.MIDI_NUM_TEMPERATURE:
                        temperature = data2 * \
                            (controls.max_temp - controls.min_temp) / \
                            127 + controls.min_temp
                        controls.temperature = limit_range(
                            temperature, controls.min_temp, controls.max_temp)
                    elif status == controls.MIDI_NUM_VOLUME:
                        volume = data2 * \
                            (controls.max_vol - controls.min_vol) / \
                            127 + controls.min_vol
                        controls.volume = limit_range(
                            volume, controls.min_vol, controls.max_vol)
                        controls.ensemble = 'NVT'
                    elif status == controls.MIDI_NUM_PRESSURE:
                        pressure = pressure_from_midi_val(
                            data2, controls.min_press, controls.max_press)
                        controls.pressure = limit_range(
                            pressure, controls.min_press, controls.max_press)
                        controls.ensemble = 'NPT'
                    elif status == controls.MIDI_NUM_NUMBEROFPARTICLES:
                        npart = int(data2 * controls.max_n / 127)
                        controls.number_of_particles = limit_range(
                            npart, controls.min_n, controls.max_n)
                    elif status == controls.MIDI_ROTATE:
                        if data2 < 65:
                            # rotate clockwise
                            mayavi_rotation_angle += mayavi_rotation_angle_step * \
                                data2
                        elif data2 >= 65:
                            # rotate counterclockwise
                            mayavi_rotation_angle -= mayavi_rotation_angle_step * \
                                (data2 - 64)
                    elif status == controls.MIDI_ZOOM:
                        if data1 == 99 and data2 == 127:
                            # zoom in
                            mayavi_zoom -= mayavi_zoom_step
                        elif data1 == 98 and data2 == 127:
                            # zoom out
                            mayavi_zoom += mayavi_zoom_step
                    # else:
                    # print("Unknown Status {0} with data1={1} and
                    # data2={2}".format(status, data1, data2))

        except Exception as e:
            print(e)

        time.sleep(0.01)


last_plotted = 0


def rotate_scene():
    global mayavi_rotation_angle

    if use_mayavi and mayavi_rotation_angle:
        # mlab.yaw(mayavi_rotation_angle)
        if mayavi_autozoom:
            mlab.view(azimuth=mayavi_rotation_angle, distance='auto')
        else:
            current_view_vals = mlab.view()
            mlab.view(azimuth=mayavi_rotation_angle,
                      elevation=current_view_vals[1],
                      distance=current_view_vals[2],
                      focalpoint=current_view_vals[3])
    mayavi_rotation_angle %= 360.


def zoom_scene():
    global mayavi_zoom, mayavi_zoom_old

    if use_mayavi:
        mlab.view(distance=mayavi_zoom)
    elif use_opengl:
        if mayavi_zoom_old < mayavi_zoom:
            vis.camera.move_backward()
            mayavi_zoom_old = mayavi_zoom
        elif mayavi_zoom_old > mayavi_zoom:
            vis.camera.move_forward()
            help(vis.camera.move_forward)
            mayavi_zoom_old = mayavi_zoom


def update_plot():
    global last_plotted

    # rotate_scene()
    zoom_scene()

    data_len = np.array([len(plt1_x_data), len(plt1_y_data),
                         len(plt2_x_data), len(plt2_y_data)]).min()
    plot1.set_xdata(plt1_x_data[:data_len])
    plot1.set_ydata(plt1_y_data[:data_len])
    plot2.set_xdata(plt2_x_data[:data_len])
    plot2.set_ydata(plt2_y_data[:data_len])

    cursor.set_offsets([plt1_y_data[data_len - 1], plt2_y_data[data_len - 1]])
    # cursor2.set_offsets([controls.temperature, controls.pressure])

    current_time = plot1.get_xdata()[-1]
    if last_plotted == current_time:
        return
    last_plotted = current_time
    plot1.axes.set_xlim(plot1.get_xdata()[0], plot1.get_xdata()[-1])
    plot1.axes.set_ylim(0.8 * plot1.get_ydata().min(),
                        1.2 * plot1.get_ydata().max())
    plot2.axes.set_xlim(plot2.get_xdata()[0], plot2.get_xdata()[-1])
    plot2.axes.set_ylim(0.8 * plot2.get_ydata().min(),
                        1.2 * plot2.get_ydata().max())
    pyplot.draw()


t = Thread(target=main_thread)
t.daemon = True
vis.register_callback(update_plot, interval=1000)
controls = Controls()
t.start()
if controls.midi_input is not None:
    t2 = Thread(target=midi_thread)
    t2.daemon = True
    t2.start()
if use_opengl:
    gui = GUI()
    vis.register_callback(gui.process_events, interval=1000)
vis.start()
