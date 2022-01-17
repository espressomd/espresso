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
Client class containing the methods used for communicating and manipulating
the system MessageID class for communication identifiers.
"""

from enum import unique, Enum, auto
from multiprocessing import Process
import os
import numpy as np


@unique
class MessageID(Enum):
    # Message identifiers
    ENERGY = auto()
    END = auto()
    MOVE_PART = auto()
    MOVE_PART_REVERT = auto()
    CHANGE_VOLUME = auto()
    CHANGE_VOLUME_REVERT = auto()
    PART_ADD = auto()
    PART_ADD_REVERT = auto()
    PART_REMOVE = auto()
    PART_REMOVE_REVERT = auto()
    INFO = auto()


class Client(Process):
    def __init__(self,
                 box_name,
                 pipe,
                 seed,
                 box_length,
                 num_part):
        super().__init__(name=box_name)
        self.init_box_length = box_length
        self.init_num_part = num_part
        self.system = None
        self.old_box_length = 0.0
        # Pipe for communicating
        self.pipe = pipe
        self.seed = seed

    @property
    def density(self):
        return len(self.system.part.all().id) / \
            float(np.prod(self.system.box_l))

    @property
    def energy(self):
        return self.system.analysis.energy()["total"]

    def send_energy(self):
        self.pipe.send([MessageID.ENERGY, self.energy])

    def end_simulation(self):
        # Additional simulation checkpointing etc.
        os._exit(os.EX_OK)

    def send_info(self):
        self.pipe.send([MessageID.INFO, self.system.box_l[0],
                        len(self.system.part.all().id)])

    def move_particle(self, DX_MAX):
        """ Move random particle inside the box """
        random_particle_id = np.random.choice(self.system.part.all().id)
        self.old_particle = self.system.part.by_id(random_particle_id)
        self.old_pos = self.old_particle.pos
        self.old_particle.pos = self.old_pos + \
            (0.5 - np.random.random(size=3)) * DX_MAX
        self.send_energy()

    def revert_move_particle(self):
        """ Revert particle move """
        self.old_particle.pos = self.old_pos
        self.old_particle = None

    def add_particle(self):
        """ Add a particle at random inside box """
        self.last_added_particle = self.system.part.add(
            pos=np.random.random(3) * self.system.box_l)
        self.send_energy()

    def revert_add_particle(self):
        """ Revert last particle add """
        self.last_added_particle.remove()
        self.last_added_particle = None

    def remove_particle(self):
        """ Remove a random particle """
        random_particle_id = np.random.choice(self.system.part.all().id)
        self.old_particle = self.system.part.by_id(
            random_particle_id).to_dict()
        self.system.part.by_id(self.old_particle["id"]).remove()
        self.send_energy()

    def revert_remove_particle(self):
        """ Revert last particle remove """
        self.system.part.add(self.old_particle)
        self.old_particle = None

    def change_volume(self, new_box_l):
        """ Change volume to new box length """
        self.old_box_length = self.system.box_l[0]
        self.system.change_volume_and_rescale_particles(new_box_l)
        self.send_energy()

    def revert_change_volume(self):
        """ Revert volume change """
        self.system.change_volume_and_rescale_particles(self.old_box_length)

    def run(self):
        """ Start while loop, wait for tasks and handle them until MessageID.END is sent """
        # Set random seed
        np.random.seed(self.seed)

        # Start the system (needs to be done here because
        # ESPResSo can only handle one system per process)
        self.init_system()

        # Send initial energy
        self.send_energy()

        messagehandler = {
            MessageID.MOVE_PART: self.move_particle,
            MessageID.MOVE_PART_REVERT: self.revert_move_particle,
            MessageID.CHANGE_VOLUME: self.change_volume,
            MessageID.CHANGE_VOLUME_REVERT: self.revert_change_volume,
            MessageID.PART_ADD: self.add_particle,
            MessageID.PART_ADD_REVERT: self.revert_add_particle,
            MessageID.PART_REMOVE: self.remove_particle,
            MessageID.PART_REMOVE_REVERT: self.revert_remove_particle,
            MessageID.INFO: self.send_info,
            MessageID.END: self.end_simulation
        }

        # callback loop, we can simply terminate the process
        # via process.terminate() in the parent process
        while True:
            # Catch errors during handling
            msg = self.pipe.recv()

            # Handle message
            messagehandler[msg[0]](*msg[1:])
