#
# Copyright (C) 2013-2021 The ESPResSo project
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
energy calculation of the systems, two instances of the Client class in
:file:`gibbs_ensemble.py` initialized (and possibly updated) in
:file:`gibbs_client_system.py` are started which each run a different instance
of ESPResSo.

The Gibbs-ensemble implemented in these scripts closely refers to chapter 8 of
Daan Frenkel and Berend Smit 'Understanding Molecular Simulation, Second edition'.
All equation references noted in this script can be found there.

Due to the cutoff and shifting of the LJ-potential the simulated points in the
phase diagram differ from the long range uncut LJ-potential. The phase diagram
of the used potential can be found in 'B. Smit, J. Chem. Phys. 96 (11), 1992,
Phase diagrams of Lennard-Jones fluids'.
"""

import multiprocessing
import numpy as np
import tqdm
import gzip
import pickle
import argparse
import logging
import enum

from gibbs_ensemble import MessageID
from client_system import Gibbs_Client


parser = argparse.ArgumentParser()
parser.add_argument('T', type=float, help="Temperature of the system")
parser.add_argument(
    '-N',
    type=int,
    default=256,
    help="Number of particles, default 256")
parser.add_argument(
    '--density',
    type=float,
    default=0.3,
    help="Init density, default 0.3")
parser.add_argument(
    '--steps',
    '-s',
    type=int,
    default=100000,
    help="Number of steps, default is 1e5")
parser.add_argument(
    '--warmup',
    '-w',
    type=int,
    default=5000,
    help="Warmup steps, default 5000")
parser.add_argument(
    '--seed',
    type=int,
    default=17,
    help="Random seed")
parser.add_argument(
    '--log',
    action='store_true',
    help="Pipe information to log file instead of stdout")
parser.add_argument(
    '--debug',
    action='store_true',
    help="Set for debugging info (use with small step number)")
args = parser.parse_args()


@enum.unique
class Moves(enum.Enum):
    # Particle moves
    MOVE_PARTICLE = 1
    EXCHANGE_VOLUME = 2
    EXCHANGE_PARTICLE = 3


# Number of client processes
NUM_OF_CLIENTS = 2

# Global number of particles
GLOBAL_NUM_PART = args.N
assert (GLOBAL_NUM_PART % 2) == 0, "The number of particles must be even"

# Init particle density
init_density = args.density

# Global volume
GLOBAL_VOLUME = GLOBAL_NUM_PART / init_density

# Init box _length
init_box_length = np.cbrt(GLOBAL_VOLUME / NUM_OF_CLIENTS)

# Temperature
kT = args.T

# Displacement
DX_MAX = 0.4
DV_MAX = 0.2

# Number of steps, warmup steps
steps = args.steps
warmup = args.warmup
# Perform about 100 system checks during the simulation
check = int(steps / 1000)

# Trial move probabilities
MOVE_CHANCE = 0.5
VOLUME_CHANCE = 0.01
EXCHANGE_CHANCE = 1 - MOVE_CHANCE - VOLUME_CHANCE

# Random seed
seed = args.seed
np.random.seed(seed)
box_seeds = np.random.randint(0, np.iinfo(np.int32).max, size=NUM_OF_CLIENTS)

# Filename
filename = f"temp_{kT}_seed_{seed}"

# Set the logging level
logging.basicConfig(
    filename=(filename + ".log" if args.log else ""),
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
    level=logging.DEBUG if args.debug else logging.INFO
)

# Number of moves, number of accepted moves
num_moves = {
    Moves.MOVE_PARTICLE: 0,
    Moves.EXCHANGE_VOLUME: 0,
    Moves.EXCHANGE_PARTICLE: 0}

num_accepted_moves = {
    Moves.MOVE_PARTICLE: 0,
    Moves.EXCHANGE_VOLUME: 0,
    Moves.EXCHANGE_PARTICLE: 0}


class Box:
    def __init__(self, box_name, pipe, seed, box_length, num_part):
        self.box_name = box_name
        self.pipe = pipe[0]
        self.energy = 0.0
        self.old_energy = 1e10
        self.box_length = box_length
        self.old_box_length = self.box_length
        self.num_part = num_part
        self.client = Gibbs_Client(
            box_name,
            pipe[1],
            seed,
            box_length,
            num_part)

    @property
    def volume(self):
        return self.box_length**3

    @property
    def density(self):
        return self.num_part / self.volume

    @property
    def diff_energy(self):
        return np.clip((self.energy - self.old_energy), -100 * kT, 100 * kT)

    def send_data(self, data):
        assert isinstance(data, list)
        logging.debug(f"Send to {self.box_name}: {data}")
        self.pipe.send(data)

    def recv_data(self):
        msg = self.pipe.recv()
        logging.debug(f"Receive from {self.box_name}: {msg}")
        return msg

    def recv_energy(self):
        msg = self.recv_data()
        # Debugging
        if msg[0] != MessageID.ENERGY:
            raise ConnectionError(
                f"Error during energy return of {self.box_name}, "
                f"got this instead: \n{msg}")
        self.old_energy = self.energy
        self.energy = msg[1]


def validate_info(boxes):
    """ Check if info returned from boxes is the same as on the parent side """
    [box.send_data([MessageID.INFO]) for box in boxes]
    for box in boxes:
        msg = box.pipe.recv()
        if msg[0] != MessageID.INFO:
            raise ConnectionError(
                f"Connection to {box.box_name} seems to be broken.")
        np.testing.assert_equal(
            box.box_length,
            msg[1],
            err_msg="Server side box length (actual) differs from client side (desired).")
        np.testing.assert_equal(
            box.num_part,
            msg[2],
            err_msg="Server side num part (actual) differs from client side (desired)")
        logging.debug(
            f"Validation correct. Values of {box.box_name}:\n"
            f"Box length:\t{box.box_length}\nNum Part:\t{box.num_part}.")


def choose_random_box(boxes):
    """ Choose a box at random (Assumes 2 Clients) """
    random_int = np.random.randint(2)
    if boxes[random_int].num_part > 0:
        return random_int
    else:
        return (random_int + 1) % 2


def perform_move_particle(boxes):
    """ Perform a particle move, check and revert if necessary """

    # Choose random box
    box = boxes[choose_random_box(boxes)]

    # Send message to Client
    box.send_data([MessageID.MOVE_PART, DX_MAX])
    box.recv_energy()

    # Debug
    logging.debug(f"Performing particle move in {box.box_name}.")

    # ---- Check move ----

    # Check move probability (see (8.3.1) in Frenkel, Smit)
    acc = np.exp(- 1. / kT * box.diff_energy)

    # Draw random number
    rand = np.random.random()

    # Check if move was valid
    if rand > acc:
        # Reject move, restore old configuration
        box.send_data([MessageID.MOVE_PART_REVERT])
        box.energy = box.old_energy
        return False

    # Debug
    logging.debug("Particle move accepted.")

    # Only if rand < acc
    return True


def perform_volume_change(boxes):
    """ Perform a random move in log(V_1/V_2) """

    # Store old box length, old volumes
    for box in boxes:
        box.old_box_length = box.box_length
    old_volume_1 = boxes[0].volume
    old_volume_2 = boxes[1].volume

    # Debugging purposes
    if args.debug:
        np.testing.assert_almost_equal(
            GLOBAL_VOLUME, old_volume_1 + old_volume_2)

    # Random move in log(V_1/V_2) (See Algorithm 18 in Frenkel, Smit)
    lnvn = np.log(old_volume_1 / old_volume_2) + \
        (np.random.random() - 0.5) * DV_MAX

    # Calculate new box lengths
    new_volume_1 = GLOBAL_VOLUME * np.exp(lnvn) / (1 + np.exp(lnvn))
    new_volume_2 = GLOBAL_VOLUME - new_volume_1
    boxes[0].box_length = np.cbrt(new_volume_1)
    boxes[1].box_length = np.cbrt(new_volume_2)

    # Debug
    logging.debug("Perform volume change:"
                  f" {boxes[0].box_name}: {old_volume_1} -> {new_volume_1};"
                  f" {boxes[1].box_name}: {old_volume_2} -> {new_volume_2}.")

    # Send new box lengths, recv new energies
    [box.send_data([MessageID.CHANGE_VOLUME, box.box_length]) for box in boxes]
    [box.recv_energy() for box in boxes]

    # ---- Check move ----

    # Calculate acceptance probability (See (8.3.3) in Frenkel, Smit)
    acc = np.power(new_volume_1 / old_volume_1, boxes[0].num_part + 1) * \
        np.power(new_volume_2 / old_volume_2, boxes[1].num_part + 1) * \
        np.exp(- 1. / kT * (boxes[0].diff_energy + boxes[1].diff_energy))

    # Draw random number
    rand = np.random.random()

    # Check if move was valid
    if rand > acc:
        # Reject move, restore old configuration
        for box in boxes:
            box.send_data([MessageID.CHANGE_VOLUME_REVERT])
            box.box_length = box.old_box_length
            box.energy = box.old_energy
        return False

    # Debug
    logging.debug("Volume change accepted.")

    # Only if rand < acc:
    return True


def perform_particle_exchange(boxes):
    """ Remove a particle of a box chosen via choose_random_box() and add it to the other """

    # Choose donor_box
    box_num = choose_random_box(boxes)
    donor_box = boxes[box_num]
    recipient_box = boxes[(box_num + 1) % 2]

    # Send moves, update num_part
    donor_box.send_data([MessageID.PART_REMOVE])
    recipient_box.send_data([MessageID.PART_ADD])
    [box.recv_energy() for box in boxes]

    # Debugging purposes
    if logging.debug:
        np.testing.assert_equal(
            donor_box.num_part +
            recipient_box.num_part,
            GLOBAL_NUM_PART)

    logging.debug(
        f"Exchange particle of {donor_box.box_name} to {recipient_box.box_name}.")

    # ---- Check move ----

    # This is the acceptance probability if the donor_box is chosen with equal
    # probability
    acc = np.exp(
        np.log(donor_box.num_part * recipient_box.volume /
               ((recipient_box.num_part + 1) * donor_box.volume))
        - 1. / kT * donor_box.diff_energy
        - 1. / kT * recipient_box.diff_energy)

    # Draw random number
    rand = np.random.random()

    # Check if move was valid
    if rand > acc:
        # Reject move, restore old configuration
        donor_box.send_data([MessageID.PART_REMOVE_REVERT])
        recipient_box.send_data([MessageID.PART_ADD_REVERT])
        donor_box.energy = donor_box.old_energy
        recipient_box.energy = recipient_box.old_energy
        return False

    donor_box.num_part -= 1
    recipient_box.num_part += 1

    # Debug
    logging.debug("Particle exchange accepted.")

    return True


def choose_move():
    """ Choose one of the trial moves with the probabilities stated at the beginning of the script """

    rand = np.random.random()

    # Choose move
    if rand < MOVE_CHANCE:
        current_move = Moves.MOVE_PARTICLE
    elif rand < MOVE_CHANCE + VOLUME_CHANCE:
        current_move = Moves.EXCHANGE_VOLUME
    else:
        current_move = Moves.EXCHANGE_PARTICLE

    return current_move


# -------------------- Start of program -------------------- #

# Lists for boxes and pipes
boxes = []
pipes = []

# Start processes
for i in range(NUM_OF_CLIENTS):
    pipes.append(multiprocessing.Pipe())
    boxes.append(Box(f"Box0{i}",
                     pipes[i],
                     box_seeds[i],
                     init_box_length,
                     int(GLOBAL_NUM_PART / 2)))

    # Start the client
    boxes[i].client.start()

    # Receive initial energy
    boxes[i].recv_energy()

logging.info("-------------------- Start of program --------------------")

logging.info(f"Simulation parameters:\nRandom seed: {args.seed},\tWarmup "
             f"steps: {warmup},\tSteps: {steps},\tFilename: {filename}")

logging.info("Warming up.")

# Warmup
for _ in (tqdm.tqdm(range(warmup)) if not args.log else range(warmup)):
    acceptance_flag = perform_move_particle(boxes)
    if acceptance_flag:
        num_accepted_moves[Moves.MOVE_PARTICLE] += 1

logging.info(
    "Particle move acceptance rate during warmup: {:.2f}%".format(
        num_accepted_moves[Moves.MOVE_PARTICLE] / warmup * 100))

# Reset the counter for particle moves
num_accepted_moves[Moves.MOVE_PARTICLE] = 0

# List of observables we want to measure during the simulation
densities_box01 = []
densities_box02 = []

logging.info("Starting simulation.")

# Run the simulation
for count in (tqdm.tqdm(range(steps)) if not args.log else range(steps)):

    current_move = choose_move()

    if current_move == Moves.MOVE_PARTICLE:
        accepted = perform_move_particle(boxes)
    elif current_move == Moves.EXCHANGE_VOLUME:
        accepted = perform_volume_change(boxes)
    else:
        accepted = perform_particle_exchange(boxes)

    # Add current move to counter
    num_moves[current_move] += 1
    if accepted:
        num_accepted_moves[current_move] += 1

    # Append simulation observables
    densities_box01.append(boxes[0].density)
    densities_box02.append(boxes[1].density)

    if count % check == 0:
        validate_info(boxes)


logging.info("Percentage of moves done:\n \
    Particle moves, Volume changes, Particle exchanges:\n \
    {:.2f}%, {:.2f}%, {:.2f}% \n \
    Acceptance rates for each:\n \
    {:.2f}%, {:.2f}%, {:.2f}%".format(
    num_moves[Moves.MOVE_PARTICLE] / steps * 100, num_moves[Moves.EXCHANGE_VOLUME] /
    steps * 100, num_moves[Moves.EXCHANGE_PARTICLE] / steps * 100,
    num_accepted_moves[Moves.MOVE_PARTICLE] / num_moves[Moves.MOVE_PARTICLE] * 100, num_accepted_moves[Moves.EXCHANGE_VOLUME] /
    num_moves[Moves.EXCHANGE_VOLUME] * 100, num_accepted_moves[Moves.EXCHANGE_PARTICLE] /
    num_moves[Moves.EXCHANGE_PARTICLE] * 100
))

logging.info("Sending close signal to boxes.")
[box.send_data([MessageID.END]) for box in boxes]

logging.info("Saving data.")
with gzip.open(filename + ".dat.gz", 'wb') as f:
    pickle.dump({
        'args': args,
        'steps': steps,
        'temperature': kT,
        'densities_box01': [float(i) for i in densities_box01],
        'densities_box02': [float(i) for i in densities_box02],
    }, f)

logging.info("-------------------- End of program --------------------")
