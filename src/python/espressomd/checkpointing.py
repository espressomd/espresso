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
from __future__ import print_function, absolute_import

from collections import OrderedDict
import sys, inspect, os, re, signal
from espressomd.utils import is_valid_type

try:
    import cPickle as pickle
except ImportError:
    import pickle


# Convenient Checkpointing for ESPResSo
class Checkpoint(object):
    """Checkpoint handling (reading and writing).

    Parameters
    ----------
    checkpoint_id : :obj:`str`
        A string identifying a specific checkpoint.
    checkpoint_path : :obj:`str`, optional
        Path for reading and writing the checkpoint.
        If not given, the CWD is used.

    """
    def __init__(self, checkpoint_id=None, checkpoint_path="."):
        # check if checkpoint_id is valid (only allow a-z A-Z 0-9 _ -)
        if not isinstance(checkpoint_id, str) or bool(re.compile(r"[^a-zA-Z0-9_\-]").search(checkpoint_id)):
            raise ValueError("Invalid checkpoint id.")

        if not isinstance(checkpoint_path, str):
            raise ValueError("Invalid checkpoint path.")

        self.checkpoint_objects = []
        self. checkpoint_signals = []
        frm = inspect.stack()[1]
        self.calling_module = inspect.getmodule(frm[0])

        checkpoint_path = os.path.join(checkpoint_path, checkpoint_id)
        self.checkpoint_dir = os.path.realpath(checkpoint_path)

        if not os.path.isdir(self.checkpoint_dir):
            os.makedirs(self.checkpoint_dir)

        # update checkpoint counter
        self.counter = 0
        while os.path.isfile(os.path.join(self.checkpoint_dir, "{}.checkpoint".format(self.counter))):
            self.counter += 1

        # init signals
        for signum in self.read_signals():
            self.register_signal(signum)


    def __getattr_submodule(self, obj, name, default):
        """
        Generalization of getattr(). __getattr_submodule(object,
        "name1.sub1.sub2", None) will return attribute sub2 if available
        otherwise None.
        
        """
        names = name.split('.')

        for i in range(len(names)-1):
            obj = getattr(obj, names[i], default)

        return getattr(obj, names[-1], default)


    def __setattr_submodule(self, obj, name, value):
        """
        Generalization of setattr(). __setattr_submodule(object,
        "name1.sub1.sub2", value) will set attribute sub2 to value. Will raise
        exception if parent modules do not exist.
        
        """
        names = name.split('.')
        tmp_obj = obj
        for i in range(len(names)-1):
            obj = getattr(obj, names[i], None)

        if obj == None:
            raise Exception("Cannot set attribute of non existing submodules: {}\nCheck the order you registered objects for checkpointing.".format(name))
        setattr(obj, names[-1], value)


    def __hasattr_submodule(self, obj, name):
        """
        Generalization of hasattr(). __hasattr_submodule(object,
        "name1.sub1.sub2") will return True if submodule sub1 has the attribute
        sub2.
        
        """
        names = name.split('.')
        for i in range(len(names)-1):
            obj = getattr(obj, names[i], None)

        return hasattr(obj, names[-1])


    def register(self, *args):
        """Register python objects for checkpointing.

        Parameters
        ----------
        args : list of :obj:`str`
            Names of python objects to be registered for checkpointing.
        
        """
        for a in args:
            if not isinstance(a, str):
                raise ValueError("The object that should be checkpointed is identified with its name given as a string.")

            #if not a in dir(self.calling_module):
            if not self.__hasattr_submodule(self.calling_module, a):
                raise KeyError("The given object '{}' was not found in the current scope.".format(a))

            if a in self.checkpoint_objects:
                raise KeyError("The given object '{}' is already registered for checkpointing.".format(a))

            self.checkpoint_objects.append(a)


    def unregister(self, *args):
        """Unregister python objects for checkpointing.

        Parameters
        ----------
        args : list of :obj:`str`
            Names of python objects to be unregistered for checkpointing.
        
        """
        for a in args:
            if not isinstance(a, str) or not a in self.checkpoint_objects:
                raise KeyError("The given object '{}' was not registered for checkpointing yet.".format(a))

            self.checkpoint_objects.remove(a)


    def get_registered_objects(self):
        """
        Returns a list of all object names that are registered for
        checkpointing.
        
        """
        return self.checkpoint_objects


    def has_checkpoints(self):
        """Check for checkpoints.

        Returns
        -------
        bool
            True if any checkpoints exist that match checkpoint_id and checkpoint_path otherwise False.
        
        """
        return self.counter > 0


    def get_last_checkpoint_index(self):
        """
        Returns the last index of the given checkpoint id. Will raise exception
        if no checkpoints are found.
        
        """
        if not self.has_checkpoints():
            raise Exception("No checkpoints found. Cannot return index for last checkpoint.")
        return self.counter-1


    def save(self, checkpoint_index=None):
        """
        Saves all registered python objects in the given checkpoint directory
        using cPickle.
        
        """
        #get attributes of registered objects
        checkpoint_data = OrderedDict()
        for obj_name in self.checkpoint_objects:
            checkpoint_data[obj_name] = self.__getattr_submodule(self.calling_module, obj_name, None)

        if checkpoint_index is None:
            checkpoint_index = self.counter
        filename = os.path.join(self.checkpoint_dir, "{}.checkpoint".format(checkpoint_index))

        tmpname = filename + ".__tmp__"
        with open(tmpname,"wb") as checkpoint_file:
            pickle.dump(checkpoint_data, checkpoint_file, -1)
        os.rename(tmpname, filename)


    def load(self, checkpoint_index=None):
        """
        Loads the python objects using (c)Pickle and sets them in the calling
        module.

        Parameters
        ----------
        checkpoint_index : :obj:`int`, optional
            If not given, the latest checkpoint_index will be used.
        
        """
        if checkpoint_index == None:
            checkpoint_index = self.get_last_checkpoint_index()

        filename = os.path.join(self.checkpoint_dir, "{}.checkpoint".format(checkpoint_index))
        with open(filename, "rb") as f:
            checkpoint_data = pickle.load(f)

        for key in checkpoint_data:
            self.__setattr_submodule(self.calling_module, key, checkpoint_data[key])
            self.checkpoint_objects.append(key)


    def __signal_handler(self, signum, frame):
        """
        Will be called when a registered signal was sent.
        
        """
        self.save()
        exit(signum)


    def read_signals(self):
        """
        Reads all registered signals from the signal file and returns a list of
        integers.
        
        """
        if not os.path.isfile(os.path.join(self.checkpoint_dir, "signals")):
            return []

        with open(os.path.join(self.checkpoint_dir, "signals"), "r") as signal_file:
            signals = signal_file.readline().strip().split()
            signals = [int(i) for i in signals] #will raise exception if signal file contains invalid entries
        return signals


    def __write_signal(self, signum=None):
        """Writes the given signal integer signum to the signal file.
        
        """
        signum = int(signum)
        if not is_valid_type(signum, int):
            raise ValueError("Signal must be an integer number.")

        signals = self.read_signals()

        if not signum in signals:
            signals.append(signum)
            signals = " ".join(str(i) for i in signals)
            with open(os.path.join(self.checkpoint_dir, "signals"), "w") as signal_file:
                signal_file.write(signals)


    def register_signal(self, signum=None):
        """Register a signal that will trigger the signal handler.

        Parameters
        ----------
        signum : :obj:`int`
            Signal to be registered.
        
        """
        if not is_valid_type(signum, int):
            raise ValueError("Signal must be an integer number.")

        if signum in self.checkpoint_signals:
            raise KeyError("The signal {} is already registered for checkpointing.".format(signum))

        signal.signal(signum, self.__signal_handler)
        self.checkpoint_signals.append(signum)
        self.__write_signal(signum)
