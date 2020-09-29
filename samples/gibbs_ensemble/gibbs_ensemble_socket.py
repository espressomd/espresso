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

import argparse
import enum
import logging as log
import matplotlib.pyplot as plt
import numpy as np
import pickle
import random
import socket
import struct
import subprocess
import time


import gibbs

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seed', type=int, nargs=1)
parser.add_argument('-S', '--steps', type=int, nargs=1)
parser.add_argument('-w', '--warm_up_steps', type=int, nargs=1)
parser.add_argument('-E', '--espresso-executable', nargs=1)
parser.add_argument('-C', '--client-script', nargs=1)
parser.add_argument('-n', '--number-of-particles', type=int, nargs=1)
parser.add_argument('-l', '--box-length', type=float, nargs=1)
parser.add_argument('-T', '--temperature', type=float, nargs=1, required=True)

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
DV_MAX = 0.1

# Monte-Carlo parameters
INIT_MOVE_CHANCE = 0.16
EXCHANGE_CHANCE = 0.5
VOLUME_CHANCE = 1.0 - INIT_MOVE_CHANCE - EXCHANGE_CHANCE

# socket parameters
HOST = 'localhost'
PORT = 31415
NUMBER_OF_CLIENTS = 2

REPORT_INTERVAL = 50  # moves


@enum.unique
class Moves(enum.Enum):
    MOVE_PARTICLE = 1
    EXCHANGE_VOLUME = 2
    EXCHANGE_PARTICLE = 3


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
    energy_old = 1.0e100
    conn = 0
    adr = 0

    def recv_energy(self):
        '''Received the energy data from the client.'''
        msg = self.recv_data()
        if msg[0] == gibbs.MessageId.ENERGY:
            self.energy = msg[1]
            return 0
        else:
            raise RuntimeError("ERROR during energy recv")

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


def select_box(n_part_0, n_part_1):
    '''
    Returns a random box index 0 or 1, where the probability is proportional
    to the number of particles in the respective boxes
    '''
    n_total = n_part_0 + n_part_1
    if random.randint(0, n_total - 1) < n_part_0:
        return 0
    else:
        return 1


def move_particle(box):
    '''
    Tries a displacement move and stores the new energy in the corresponding box
    '''
    box.send_data(pickle.dumps([gibbs.MessageId.MOVE_PART]))
    box.recv_energy()


def exchange_volume(boxes, global_volume):
    '''
    Tries a volume exchange move and stores the new energy in the boxes
    '''

    rand_box = random.randint(0, 1)
    if boxes[rand_box].n_particles == 0:
        rand_box = (rand_box + 1) % 2

    lead_box = boxes[rand_box]
    other_box = boxes[(rand_box + 1) % 2]

    lead_box.box_l_old = lead_box.box_l
    other_box.box_l_old = other_box.box_l

    vol2 = global_volume - lead_box.box_l**3
    lnvn = np.log(lead_box.box_l**3 / vol2) + \
        (random.random() - 0.5) * DV_MAX
    vol1 = global_volume * np.exp(lnvn) / (1 + np.exp(lnvn))
    vol2 = global_volume - vol1
    log.info("exchange volume %g->%g", lead_box.box_l, np.cbrt(vol1))

    lead_box.box_l = np.cbrt(vol1)
    lead_box.send_data(pickle.dumps(
        [gibbs.MessageId.CHANGE_VOLUME, lead_box.box_l]))

    other_box.box_l = np.cbrt(vol2)
    other_box.send_data(pickle.dumps(
        [gibbs.MessageId.CHANGE_VOLUME, other_box.box_l]))

    lead_box.recv_energy()
    other_box.recv_energy()


def exchange_particle(boxes):
    '''
    Tries a particle exchange move and stores the new energy in the boxes
    '''
    rand_box = random.randint(0, 1)
    if boxes[rand_box].n_particles == 0:
        rand_box = (rand_box + 1) % 2

    log.info("exchange particle %d", rand_box)
    source_box = boxes[rand_box]
    dest_box = boxes[(rand_box + 1) % 2]

    source_box.n_particles -= 1
    dest_box.n_particles += 1

    source_box.send_data(pickle.dumps(
        [gibbs.MessageId.EXCHANGE_PART_REMOVE]))
    dest_box.send_data(pickle.dumps([gibbs.MessageId.EXCHANGE_PART_ADD, 1]))

    source_box.recv_energy()
    dest_box.recv_energy()

    return source_box, dest_box


def check_make_move(box, inner_potential):
    '''
    Returns whether the last displacement move was valid or not and reverts the changes
    if it was invalid.
    '''
    if random.random() > inner_potential:
        log.info("revert move")
        box.send_data(pickle.dumps([gibbs.MessageId.MOVE_PART_REVERT]))
        box.energy = box.energy_old
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
        log.info("revert exchange volume")
        boxes[0].send_data(pickle.dumps(
            [gibbs.MessageId.CHANGE_VOLUME_REVERT, boxes[0].box_l_old]))
        boxes[0].box_l = boxes[0].box_l_old
        boxes[1].send_data(pickle.dumps(
            [gibbs.MessageId.CHANGE_VOLUME_REVERT, boxes[1].box_l_old]))
        boxes[1].box_l = boxes[1].box_l_old

        boxes[0].energy = boxes[0].energy_old
        boxes[1].energy = boxes[1].energy_old

        return False
    return True


def check_exchange_particle(source_box, dest_box, inner_potential):
    '''
    Returns whether the last particle exchange move was valid or not and reverts the changes
    if it was invalid.
    '''
    exchange_factor = (source_box.n_particles) / \
        (dest_box.n_particles + 1.0) * \
        dest_box.box_l**3 / source_box.box_l**3
    if random.random() > exchange_factor * inner_potential:
        log.info("revert exchange particle")
        source_box.n_particles += 1
        dest_box.n_particles -= 1
        source_box.send_data(pickle.dumps(
            [gibbs.MessageId.EXCHANGE_PART_REMOVE_REVERT]))
        dest_box.send_data(pickle.dumps(
            [gibbs.MessageId.EXCHANGE_PART_ADD_REVERT]))

        source_box.energy = source_box.energy_old
        dest_box.energy = dest_box.energy_old
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

    arguments = ["--seed", str(random.randint(0, np.iinfo(np.int32).max)), 
                 "--port", str(PORT)]
    subprocess.Popen([espresso_executable] + [client_script] +
                     arguments)

    boxes[i].conn, boxes[i].adr = s.accept()
    boxes[i].send_data(pickle.dumps([gibbs.MessageId.START]))

    boxes[i].recv_energy()

# set initial volume and particles
for box in boxes:
    box.send_data(pickle.dumps(
        [gibbs.MessageId.CHANGE_VOLUME, box.box_l]))
    box.recv_energy()
    box.send_data(pickle.dumps(
        [gibbs.MessageId.EXCHANGE_PART_ADD, box.n_particles]))
    box.recv_energy()

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


# only do displacements during the warm up
move_chance = 1.0
# Monte-Carlo loop
tick = time.time()
for i in range(steps):
    # print("\rIntegrating: {0:3.0f}%".format(
    #    100 * i / steps), end='', flush=True)

    # warm up ends after 1000 steps
    if i == warm_up_steps:
        move_chance = INIT_MOVE_CHANCE

    # rand defines which move to make.
    rand = random.random()

    # choose step and send the command to execute to the clients, then receive
    # the energies.
    if rand <= move_chance:
        current_move = Moves.MOVE_PARTICLE
        idx = select_box(boxes[0].n_particles, boxes[1].n_particles)
        log.info("move particle in box %d", idx)
        box = boxes[idx]
        move_particle(box)
        move_steps_tried += 1

    elif rand <= move_chance + VOLUME_CHANCE:
        current_move = Moves.EXCHANGE_VOLUME
        exchange_volume(boxes, global_volume)
        volume_steps_tried += 1

    else:
        current_move = Moves.EXCHANGE_PARTICLE
        source_box, dest_box = exchange_particle(boxes)
        exchange_steps_tried += 1

    # test if the move will be accepted and undo the last step if it was not
    # accepted.
    delta_energy = boxes[0].energy + boxes[1].energy - \
        boxes[0].energy_old - boxes[1].energy_old
    # limitation to delta_energy to circumvent calculating the exponential of
    # too large inner potentials, which could cause some error messages.
    if delta_energy < -10.0 * kT:
        delta_energy = -10.0 * kT

    log.info("dE: %g", delta_energy)
    inner_potential = np.exp(-1.0 / kT * delta_energy)
    if current_move == Moves.MOVE_PARTICLE:
        accepted = check_make_move(box, inner_potential)
    elif current_move == Moves.EXCHANGE_VOLUME:
        accepted = check_exchange_volume(boxes, inner_potential)
    elif current_move == Moves.EXCHANGE_PARTICLE: 
        accepted = check_exchange_particle(
            source_box, dest_box, inner_potential)
    else:
        raise Exception("Unkhandled move type " + current_move)

    # accounting
    boxes[0].energy_old = boxes[0].energy
    boxes[1].energy_old = boxes[1].energy
    if accepted:
        log.info("accepted")
        # keep the changes.
        accepted_steps += 1
        if current_move == Moves.MOVE_PARTICLE:
            accepted_steps_move += 1
        elif current_move == Moves.EXCHANGE_VOLUME:
            accepted_steps_volume += 1
        elif current_move == Moves.EXCHANGE_PARTICLE:
            accepted_steps_exchange += 1
        else: 
            raise Exception("Unhandled move type " + current_move)

        densities[0].append(boxes[0].n_particles / boxes[0].box_l**3)
        densities[1].append(boxes[1].n_particles / boxes[1].box_l**3)
        list_indices.append(i)

        if i % REPORT_INTERVAL == 0 and i > 0:
            tock = time.time()
            print("step %d, densities %.3f %.3f, %.2f sec / move" % 
                  (i, densities[0][-1], densities[1][-1], 
                   (tock - tick) / REPORT_INTERVAL))
            assert boxes[0].n_particles + \
                boxes[1].n_particles == global_num_particles
            assert abs(boxes[0].box_l**3 + boxes[1].box_l**3 - global_volume) \
                < 1E-8 * global_volume
            tick = time.time()


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
    boxes[i].send_data(pickle.dumps([gibbs.MessageId.END]))
s.close()
