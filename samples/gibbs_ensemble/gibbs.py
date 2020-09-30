import enum
import pickle
import socket
import struct

import numpy as np


@enum.unique
class MessageId(enum.Enum):
    START = 0
    END = 1
    MOVE_PART = 2
    MOVE_PART_REVERT = 21
    CHANGE_VOLUME = 3
    CHANGE_VOLUME_REVERT = 31
    EXCHANGE_PART_ADD = 4
    EXCHANGE_PART_ADD_REVERT = 41
    EXCHANGE_PART_REMOVE = 5
    EXCHANGE_PART_REMOVE_REVERT = 51
    ENERGY = 6
    CONSISTENCY_CHECK = 7


class Client:
    '''
    Base class for a Gibbs ensemble client. Each client runs a single
    system, executes requested changes and returns energies.
    '''

    def __init__(self, system):
        self._system = system

    def random_particle(self):
        id = np.random.choice(self._system.part[:].id)
        return self._system.part[id]

    def move_particle(self):
        '''Moves a particle to a new random position. The old position is saved in _old_part_position.'''
        p = self.random_particle()
        self._last_modified_pid = p.id
        self._saved_pos = p.pos
        p.pos = np.random.random(3) * self._system.box_l

    def move_particle_revert(self):
        '''Revert the last movement step.'''
        self._system.part[self._last_modified_pid].pos = self._saved_pos

    def particle_state(self, p):
        """This implementation only saves position"""
        return {'id': p.id, 'pos': p.pos}

    def remove_particle(self):
        '''Removes a random particle. The old position and index are saved'''
        p = self.random_particle()
        self._saved_particle_state = self.particle_state(p)
        p.remove()

    def remove_particle_revert(self):
        '''Revert the last remove particle step.'''
        self._system.part.add(self._saved_particle_state)

    def change_volume_and_rescale_particles(self, new_box_l):
        '''Changes the volume and rescales the particle positions.'''
        self._system.change_volume_and_rescale_particles(new_box_l)

    def insert_particle(self, n):
        '''Inserts n particles at random position.'''
        self.last_added_particles = self._system.part.add(
            pos=np.random.random((n, 3)) * self._system.box_l)

    def remove_last_added_particle(self):
        '''Removes the last added particle (reverts insert particle).'''
        for p in self.last_added_particles: p.remove()

    def energy(self):
        return self._system.analysis.energy()["total"]

    def handle_move_particle(self):
        self.move_particle()
        send_data(self._socket,
                  [MessageId.ENERGY, self.energy()])

    def handle_move_particle_revert(self):
        self.move_particle_revert()

    def handle_change_volume(self, box_l):
        self.change_volume_and_rescale_particles(box_l)
        send_data(self._socket,
                  [MessageId.ENERGY, self.energy()])

    def handle_change_volume_revert(self, box_l): 
        self.change_volume_and_rescale_particles(box_l)

    def handle_add_particle(self, n):
        self.insert_particle(n)
        send_data(self._socket,
                  [MessageId.ENERGY, self.energy()])

    def handle_add_particle_revert(self):
        self.remove_last_added_particle()

    def handle_remove_particle(self):
        self.remove_particle()
        send_data(self._socket,
                  [MessageId.ENERGY, self.energy()])

    def handle_remove_particle_revert(self):
        self.remove_particle_revert()

    def handle_consistency_check(self):
        send_data(self._socket,
                  [len(self._system.part),
                   self._system.box_l[0]])

    def run(self, host, port):
        # mapping between messages and handler functions
        message_handlers = {
            MessageId.MOVE_PART: self.handle_move_particle,
            MessageId.MOVE_PART_REVERT: self.handle_move_particle_revert,

            MessageId.CHANGE_VOLUME: self.handle_change_volume,
            MessageId.CHANGE_VOLUME_REVERT: self.handle_change_volume_revert,

            MessageId.EXCHANGE_PART_ADD: self.handle_add_particle,
            MessageId.EXCHANGE_PART_ADD_REVERT: self.handle_add_particle_revert,

            MessageId.EXCHANGE_PART_REMOVE: self.handle_remove_particle,
            MessageId.EXCHANGE_PART_REMOVE_REVERT: self.handle_remove_particle_revert,

            MessageId.CONSISTENCY_CHECK: self.handle_consistency_check}

        # init socket
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect((host, port))
        msg = recv_data(self._socket)
        if msg[0] != MessageId.START:
            raise Exception("Invalid message received from controller")

        # send the initial energy
        energy = self.energy()
        send_data(self._socket, [MessageId.ENERGY, energy])

        while msg[0] != MessageId.END:
            # receive command to execute next step
            msg = recv_data(self._socket)
            if msg[0] == MessageId.END:
                break

            # Handle remaining messages according to the registered message handlers
            # Pass any arguments to the handler.
            message_handlers[msg[0]](*msg[1:])

        # close the socket
        self._socket.close()


MAX_MSG_SIZE = 128 


def recv_data(socket):
    '''Receives data and return it.'''
    # The first 4 bytes encode the length of the messages received.
    return pickle.loads(socket.recv(MAX_MSG_SIZE))


def send_data(socket, data):
    '''Send the data packet.'''
    # Get rid of data overhead caused by numpy data types
    converted = [float(i) if isinstance(i, np.float) else i for i in data]
    pickled = pickle.dumps(converted)
    assert len(pickled) < MAX_MSG_SIZE
    socket.send(pickled)
