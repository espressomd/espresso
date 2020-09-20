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
Client part of the Gibbs ensemble simulation. This script handles the
simulation boxes and communicates the energies to the host. The Monte-Carlo
part of the simulation is done by the :file:`gibbs_ensemble_socket.py` script.
"""

import espressomd
import socket
import numpy as np
import pickle
import struct
import argparse
import gibbs


espressomd.assert_features("LENNARD_JONES")

seed = None
init_box_l = None
particle_number = None

lj_epsilon = 1.0
lj_sigma = 1.0
lj_cutoff = 2.5
lj_shift = 0

HOST = 'localhost'
PORT = 31415

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--number-particles', type=int, nargs=1)
parser.add_argument('-s', '--seed', type=int, nargs=1)
parser.add_argument('-bl', '--box-length', type=float, nargs=1)

args = parser.parse_args()

if args.number_particles:
    particle_number = args.number_particles[0]
if args.seed:
    seed = args.seed[0]
if args.box_length:
    init_box_l = args.box_length[0]


class Manipulator():
    '''
    Manipulator class, which has all the system changing functions.
    '''

    def __init__(self, system):
        self._system = system
        self._old_part_position = []
        self._old_part_idx = None

    def move_particle(self):
        '''Moves a particle to a new random position. The old position is saved in _old_part_position.'''
        sel = self._system.part.select(lambda p: p.id > -1).id_selection
        self._old_part_idx = np.random.choice(sel)
        self._old_part_position = self._system.part[self._old_part_idx].pos
        self._system.part[self._old_part_idx].pos = np.random.rand(3) * box_l

    def move_particle_revert(self):
        '''Revert the last movement step.'''
        self._system.part[self._old_part_idx].pos = self._old_part_position

    def remove_particle(self):
        '''Removes a random particle. The old position and index are saved in _old_part_position and _old_part_idx.'''
        sel = self._system.part.select(lambda p: p.id > -1).id_selection
        self._old_part_idx = np.random.choice(sel)
        self._old_part_position = self._system.part[self._old_part_idx].pos
        self._system.part[self._old_part_idx].remove()

    def remove_particle_revert(self):
        '''Revert the last remove particle step.'''
        self._system.part.add(pos=self._old_part_position)

    def change_volume_and_rescale_particles(self, new_box_l):
        '''Changes the volume and rescales the particle positions.'''
        self._system.change_volume_and_rescale_particles(new_box_l)

    def insert_particle(self):
        '''Inserts a particle at a random position.'''
        self._system.part.add(pos=np.random.rand(3) * self._system.box_l)

    def remove_last_added_particle(self):
        '''Removes the last added particle (reverts insert particle).'''
        self._system.part[self._system.part.highest_particle_id].remove()

    def energy(self):
        density = len(self._system.part) / self._system.volume()
        E_sim = self._system.analysis.energy()["total"]
        E_tail = calc_tail_correction(density)
        E_shift = calc_shift_correction(density)
        return E_sim + E_tail + E_shift


def calc_tail_correction(density):
    '''
    Calculates the tail correction to the energies of the box.
    '''
    # eq 3.2.5
    return 8.0 / 3.0 * np.pi * density * lj_epsilon * \
        lj_sigma**3 * (1.0 / 3.0 * np.power(lj_cutoff / lj_sigma, -9) -
                       np.power(lj_cutoff / lj_sigma, -3))


def calc_shift_correction(density):
    '''
    Calculates the shift correction to the energies of the box.
    '''
    # difference in the potential integrated from 0 to cutoff distance
    return -8.0 / 3.0 * np.pi * density * \
        lj_epsilon * np.power(lj_cutoff, 3) * 4.0 * lj_shift


def recv_data(socket):
    '''Receives data and return it.'''
    # The first 4 bytes encode the length of the messages received.
    buf = b''
    while len(buf) < 4:
        buf += socket.recv(4 - len(buf))
    length = struct.unpack('!I', buf)[0]
    msg = pickle.loads(socket.recv(length))
    return msg


def send_data(data, socket):
    '''Send the data packet.'''
    # The first 4 bytes encode the length of the messages sent.
    length = struct.pack('>I', len(data))
    packet = length + data
    socket.send(packet)


# set random seed
np.random.seed(seed)

# init socket
socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
socket.connect((HOST, PORT))
msg = recv_data(socket)

# init system
box_l = init_box_l
system = espressomd.System(box_l=[box_l, box_l, box_l])
system.cell_system.set_n_square()
system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=lj_epsilon,
                                                       sigma=lj_sigma,
                                                       cutoff=lj_cutoff,
                                                       shift=lj_shift)
system.time_step = 0.01
manipulator = Manipulator(system)

# places the particles
for i in range(particle_number):
    system.part.add(pos=np.random.rand(3) * box_l, type=0, id=i)

# send the initial energy
energy = manipulator.energy()
send_data(pickle.dumps([gibbs.MessageId.ENERGY, energy]), socket)

while msg[0] != gibbs.MessageId.END:
    # receive command to execute next step
    msg = recv_data(socket)
    if msg[0] == gibbs.MessageId.END:
        break
    elif msg[0] == gibbs.MessageId.MOVE_PART:
        manipulator.move_particle()
    elif msg[0] == gibbs.MessageId.CHANGE_VOLUME:
        box_l = msg[1]
        manipulator.change_volume_and_rescale_particles(box_l)
    elif msg[0] == gibbs.MessageId.EXCHANGE_PART_ADD:
        manipulator.insert_particle()
    elif msg[0] == gibbs.MessageId.EXCHANGE_PART_REMOVE:
        manipulator.remove_particle()
    elif msg[0] == gibbs.MessageId.MOVE_PART_REVERT:
        manipulator.move_particle_revert()
    elif msg[0] == gibbs.MessageId.EXCHANGE_PART_ADD_REVERT:
        manipulator.remove_last_added_particle()
    elif msg[0] == gibbs.MessageId.EXCHANGE_PART_REMOVE_REVERT:
        manipulator.remove_particle_revert()

    # calculation energy and send it to the host
    energy = manipulator.energy()
    send_data(pickle.dumps([gibbs.MessageId.ENERGY, energy]), socket)

# closing the socket
socket.close()
