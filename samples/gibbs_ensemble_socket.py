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
Simulate a Gibbs-ensemble of a Lennard-Jones fluid at a fixed temperature.
This script does the Monte-Carlo part of the Gibbs-ensemble, however for the
energy calculation of the systems, two instances of the :file:`gibbs_ensemble_client.py`
script are executed, which each run a different instance of ESPResSo.

The Gibbs-ensemble implemented in these scripts closely refers to chapter 8 of
Daan Fenkel and Berend Smit 'Understanding Molecular Simulation, Second edition'.
All equation references noted in this script can be found there.

Due to the cutoff and shifting of the LJ-potential the simulated points in the
phase diagram differ from the long range uncut LJ-potential. The phase diagram
of the used potential can be found in 'B. Smit, J. Chem. Phys. 96 (11), 1992,
Phase diagrams of Lennard-Jones fluids'.
"""

import socket
import numpy as np
import pickle
import subprocess
import struct
import random
import matplotlib.pyplot as plt
import argparse

from espressomd import assert_features
assert_features("LENNARD_JONES")

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, nargs=1)
parser.add_argument('-S', '--steps', type=int, nargs=1)
parser.add_argument('-w', '--warm_up_steps', type=int, nargs=1)
parser.add_argument('-E', '--espresso-executable', nargs=1)
parser.add_argument('-C', '--client-script', nargs=1)
parser.add_argument('-n', '--number-of-particles', type=int, nargs=1)
parser.add_argument('-l', '--box-length', type=float, nargs=1)
parser.add_argument('-T', '--temperature', type=float, nargs=1)

args = parser.parse_args()

# system parameters
seed = None
steps = 50000
warm_up_steps = 1000
# number of particles in both boxes combined
global_num_particles = 256
# starting box length
init_box_l = 7.55
# temperature
kT = 1.15
# maximum of the volume exchanged in the logarithmic space
DV_MAX = 0.5
PARTICLE_RADIUS = 0.5

# LJ-parameters
LJ_EPSILON = 1.0
LJ_SIGMA = 2.0 * PARTICLE_RADIUS
LJ_CUTOFF = 2.5
LJ_SHIFT = 0.0

# Monte-Carlo parameters
INIT_MOVE_CHANCE = 0.16
EXCHANGE_CHANCE = 0.8
VOLUME_CHANCE = 1.0 - INIT_MOVE_CHANCE - EXCHANGE_CHANCE

# socket parameters
HOST = 'localhost'
PORT = 31415
NUMBER_OF_CLIENTS = 2

# Message identifiers
MSG_START = 0
MSG_END = 1 
MSG_MOVE_PART = 2
MSG_MOVE_PART_REVERT = 21
MSG_CHANGE_VOLUME = 3
MSG_EXCHANGE_PART_ADD = 4
MSG_EXCHANGE_PART_ADD_REVERT = 41
MSG_EXCHANGE_PART_REMOVE = 5
MSG_EXCHANGE_PART_REMOVE_REVERT = 51
MSG_ENERGY = 6

# script locations
espresso_executable = "../pypresso"
client_script = "./gibbs_ensemble_client.py"

if args.seed:
    seed = args.seed[0]
if args.steps:
    steps = args.steps[0]
if args.warm_up_steps:
    warm_up_steps = args.warm_up_steps[0]
if args.espresso_executable:
    espresso_executable = args.espresso_executable[0]
if args.client_script:
    client_script = args.client_script[0]
if args.number_of_particles:
    global_num_particles = args.number_of_particles[0]
if args.temperature:
    kT = args.temperature[0]
if args.box_length:
    init_box_l = args.box_length[0]

global_volume = 2.0 * init_box_l**3


class Box:
    '''
    Box class, which contains the data of one system and controls the
    connection to the respective client.
    '''
    box_l = init_box_l
    box_l_old = init_box_l
    n_particles = int(global_num_particles / 2)
    energy = 0.0
    energy_corrected = 0.0
    energy_old = 1.0e100
    conn = 0
    adr = 0

    def recv_energy(self):
        '''Received the energy data from the client.'''
        msg = self.recv_data()
        if msg[0] == MSG_ENERGY:
            self.energy = msg[1]
            return 0
        else:
            print("ERROR during energy recv")
            return 1

    def recv_data(self):
        '''Received the data send from the client.'''
        # The first 4 bytes encode the length of the messages received.
        buf = b''
        while len(buf) < 4:
            buf += self.conn.recv(4 - len(buf))
        length = struct.unpack('!I', buf)[0]
        d = self.conn.recv(length)
        msg = pickle.loads(d)
        return(msg)

    def send_data(self, data):
        '''Send the data packet to the client.'''
        # The first 4 bytes encode the length of the messages sent.
        length = struct.pack('>I', len(data))
        packet = length + data
        self.conn.send(packet)


def calc_tail_correction(box, lj_epsilon, lj_sigma, lj_cutoff):
    '''
    Calculates the tail correction to the energies of the box.
    '''
    # eq 3.2.5
    return 8.0 / 3.0 * np.pi * box.n_particles / box.box_l**3 * lj_epsilon * \
        lj_sigma**3 * (1.0 / 3.0 * np.power(lj_cutoff / lj_sigma, -9) -
                       np.power(lj_cutoff / lj_sigma, -3))


def calc_shift_correction(box, lj_epsilon, lj_cutoff, lj_shift):
    '''
    Calculates the shift correction to the energies of the box.
    '''
    # difference in the potential integrated from 0 to cutoff distance
    return -8.0 / 3.0 * np.pi * box.n_particles / box.box_l**3 * \
        lj_epsilon * np.power(lj_cutoff, 3) * 4.0 * lj_shift


def move_particle(boxes, global_num_particles):
    '''
    Tries a displacement move and stores the new energy in the corresponding box
    '''
    if random.randint(0, global_num_particles - 1) < boxes[0].n_particles:
        rand_box = 0
    else:
        rand_box = 1
    boxes[rand_box].send_data(pickle.dumps([MSG_MOVE_PART, 0]))
    boxes[rand_box].recv_energy()
    return rand_box


def exchange_volume(boxes, global_volume):
    '''
    Tries a volume exchange move and stores the new energy in the boxes
    '''
    boxes[0].box_l_old = boxes[0].box_l
    boxes[1].box_l_old = boxes[1].box_l

    # calculate the exchanged volume
    rand_box = random.randint(0, NUMBER_OF_CLIENTS - 1)
    vol2 = global_volume - boxes[rand_box].box_l**3
    lnvn = np.log(boxes[rand_box].box_l**3 / vol2) + \
        (random.random() - 0.5) * DV_MAX
    vol1 = global_volume * np.exp(lnvn) / (1 + np.exp(lnvn))
    vol2 = global_volume - vol1
    boxes[rand_box].box_l = np.cbrt(vol1)
    boxes[rand_box].send_data(pickle.dumps(
        [MSG_CHANGE_VOLUME, boxes[rand_box].box_l]))

    boxes[(rand_box + 1) % 2].box_l = np.cbrt(vol2)
    boxes[(rand_box + 1) % 2].send_data(pickle.dumps([MSG_CHANGE_VOLUME,
                                                      boxes[(rand_box + 1) % 2].box_l]))

    boxes[rand_box].recv_energy()
    boxes[(rand_box + 1) % 2].recv_energy()
    return rand_box


def exchange_particle(boxes):
    '''
    Tries a particle exchange move and stores the new energy in the boxes
    '''
    rand_box = random.randint(0, 1)
    if boxes[rand_box].n_particles == 0:
        rand_box = (rand_box + 1) % 2
    boxes[rand_box].n_particles -= 1
    boxes[(rand_box + 1) % 2].n_particles += 1
    boxes[rand_box].send_data(pickle.dumps([MSG_EXCHANGE_PART_REMOVE, 0]))
    boxes[(rand_box + 1) %
          2].send_data(pickle.dumps([MSG_EXCHANGE_PART_ADD, 0]))

    boxes[rand_box].recv_energy()
    boxes[(rand_box + 1) % 2].recv_energy()
    return rand_box


def check_make_move(boxes, inner_potential, rand_box):
    '''
    Returns whether the last displacement move was valid or not and reverts the changes
    if it was invalid.
    '''
    if random.random() > inner_potential:
        boxes[rand_box].send_data(pickle.dumps([MSG_MOVE_PART_REVERT, 0]))
        boxes[rand_box].recv_energy()
        return False
    return True


def check_exchange_volume(boxes, inner_potential):
    '''
    Returns whether the last volume exchange move was valid or not and reverts the changes
    if it was invalid.
    '''
    volume_factor = \
        (boxes[0].box_l**3 / boxes[0].box_l_old**3)**(boxes[0].n_particles + 1) * \
        (boxes[1].box_l**3 / boxes[1].box_l_old **
         3)**(boxes[1].n_particles + 1)
    if random.random() > volume_factor * inner_potential:
        boxes[0].send_data(pickle.dumps(
            [MSG_CHANGE_VOLUME, boxes[0].box_l_old]))
        boxes[0].box_l = boxes[0].box_l_old
        boxes[1].send_data(pickle.dumps(
            [MSG_CHANGE_VOLUME, boxes[1].box_l_old]))
        boxes[1].box_l = boxes[1].box_l_old

        boxes[0].recv_energy()
        boxes[1].recv_energy()
        return False
    return True


def check_exchange_particle(boxes, inner_potential, rand_box):
    '''
    Returns whether the last particle exchange move was valid or not and reverts the changes
    if it was invalid.
    '''
    exchange_factor = (boxes[rand_box].n_particles) / \
        (boxes[(rand_box + 1) % 2].n_particles + 1.0) * \
        boxes[(rand_box + 1) % 2].box_l**3 / boxes[rand_box].box_l**3
    if random.random() > exchange_factor * inner_potential:
        boxes[rand_box].n_particles += 1
        boxes[(rand_box + 1) % 2].n_particles -= 1
        boxes[rand_box].send_data(pickle.dumps(
            [MSG_EXCHANGE_PART_REMOVE_REVERT, 0]))
        boxes[(rand_box + 1) %
              2].send_data(pickle.dumps([MSG_EXCHANGE_PART_ADD_REVERT, 0]))

        boxes[rand_box].recv_energy()
        boxes[(rand_box + 1) % 2].recv_energy()
        return False
    return True



# init socket
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
s.bind((HOST, PORT))
s.listen(NUMBER_OF_CLIENTS)

boxes = []

random.seed(seed)

# start clients and connections
for i in range(NUMBER_OF_CLIENTS):
    boxes.append(Box())
    lj_arguments = ["-lj", str(LJ_EPSILON), str(LJ_SIGMA), str(LJ_CUTOFF),
                    str(LJ_SHIFT)]
    arguments = ["-n", str(boxes[i].n_particles), "-s",
                 str(random.randint(0, np.iinfo(np.int32).max)), "-bl",
                 str(boxes[i].box_l)]
    subprocess.Popen([espresso_executable] + [client_script] +
                     arguments + lj_arguments)

    boxes[i].conn, boxes[i].adr = s.accept()
    boxes[i].send_data(pickle.dumps([MSG_START, 0]))

    boxes[i].recv_energy()
    boxes[i].energy_old = boxes[i].energy + \
        calc_tail_correction(boxes[i], LJ_EPSILON, LJ_SIGMA, LJ_CUTOFF) + \
        calc_shift_correction(boxes[i], LJ_EPSILON, LJ_CUTOFF, LJ_SHIFT)

# observables
densities = [[], []]
list_indices = []
accepted_steps = 0
accepted_steps_move = 1
accepted_steps_volume = 1
accepted_steps_exchange = 1
move_steps_tried = 1
volume_steps_tried = 1
exchange_steps_tried = 1

# rand_box defines on which box the move acts on.
rand_box = 0

# only do displacements during the warm up
move_chance = 1.0
# Monte-Carlo loop
for i in range(steps):
    print("\rIntegrating: {0:3.0f}%".format(
        100 * i / steps), end='', flush=True)

    # warm up ends after 1000 steps
    if i == warm_up_steps:
        move_chance = INIT_MOVE_CHANCE

    # rand defines which move to make.
    rand = random.random()

    # choose step and send the command to execute to the clients, then receive
    # the energies.
    if rand <= move_chance:
        rand_box = move_particle(boxes, global_num_particles)
        move_steps_tried += 1

    elif rand <= move_chance + VOLUME_CHANCE:
        rand_box = exchange_volume(boxes, global_volume)
        volume_steps_tried += 1

    else:
        rand_box = exchange_particle(boxes)
        exchange_steps_tried += 1

    # calculate the correction energies of the lj tail and shift.
    if rand <= move_chance:
        boxes[rand_box].energy_corrected = boxes[rand_box].energy + \
            calc_tail_correction(boxes[rand_box], LJ_EPSILON, LJ_SIGMA, LJ_CUTOFF) + \
            calc_shift_correction(
            boxes[rand_box],
            LJ_EPSILON,
            LJ_CUTOFF,
            LJ_SHIFT)
        boxes[(rand_box + 1) % 2].energy_corrected = \
            boxes[(rand_box + 1) % 2].energy_old

    else:
        boxes[0].energy_corrected = boxes[0].energy + \
            calc_tail_correction(boxes[0], LJ_EPSILON, LJ_SIGMA, LJ_CUTOFF) + \
            calc_shift_correction(boxes[0], LJ_EPSILON, LJ_CUTOFF, LJ_SHIFT)
        boxes[1].energy_corrected = boxes[1].energy + \
            calc_tail_correction(boxes[1], LJ_EPSILON, LJ_SIGMA, LJ_CUTOFF) + \
            calc_shift_correction(boxes[1], LJ_EPSILON, LJ_CUTOFF, LJ_SHIFT)

    # test if the move will be accepted and undo the last step if it was not
    # accepted.
    delta_energy = boxes[0].energy_corrected + boxes[1].energy_corrected - \
        boxes[0].energy_old - boxes[1].energy_old
    # limitation to delta_energy to circumvent calculating the exponential of
    # too large inner potentials, which could cause some error messages.
    if delta_energy < -10.0 * kT:
        delta_energy = -10.0 * kT

    inner_potential = np.exp(-1.0 / kT * delta_energy)
    accepted = True

    if rand <= move_chance:
        accepted = check_make_move(boxes, inner_potential, rand_box)
    elif rand <= move_chance + VOLUME_CHANCE:
        accepted = check_exchange_volume(boxes, inner_potential)
    else:
        accepted = check_exchange_particle(boxes, inner_potential, rand_box)

    if accepted:
        # keep the changes.
        boxes[0].energy_old = boxes[0].energy_corrected
        boxes[1].energy_old = boxes[1].energy_corrected
        accepted_steps += 1
        if rand <= move_chance:
            accepted_steps_move += 1
        elif rand <= move_chance + VOLUME_CHANCE:
            accepted_steps_volume += 1
        else:
            accepted_steps_exchange += 1
        densities[0].append(boxes[0].n_particles / boxes[0].box_l**3)
        densities[1].append(boxes[1].n_particles / boxes[1].box_l**3)

        list_indices.append(i)


print("Acceptance rate global:\t {}".format(accepted_steps / float(steps)))
print("Acceptance rate moving:\t {}".format(
    accepted_steps_move / float(move_steps_tried)))
print("Acceptance rate volume exchange:\t {}".format(
    accepted_steps_volume / float(volume_steps_tried)))
print("Acceptance rate particle exchange:\t {}".format(
    accepted_steps_exchange / float(exchange_steps_tried)))


plt.figure()
plt.ylabel('density')
plt.xlabel('number of steps')
plt.plot(list_indices[100:], densities[0][100:], label="box 0")
plt.plot(list_indices[100:], densities[1][100:], label="box 1")
plt.legend()
plt.show()

# closing the socket
for i in range(NUMBER_OF_CLIENTS):
    boxes[i].send_data(pickle.dumps([MSG_END, 0]))
s.close()
