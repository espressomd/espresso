#
# Copyright (C) 2013-2024 The ESPResSo project
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
import collections
import inspect
import pickle
import re
import signal
import pathlib
from . import utils
from . import script_interface


class Checkpoint:

    """Checkpoint handling (reading and writing).

    Parameters
    ----------
    checkpoint_id : :obj:`str`
        A string identifying a specific checkpoint.
    checkpoint_path : :obj:`str` or :obj:`pathlib.Path`, optional
        Path for reading and writing the checkpoint.
        If not given, the current working directory is used.

    """

    def __init__(self, checkpoint_id=None, checkpoint_path="."):
        # check if checkpoint_id is valid (only allow a-z A-Z 0-9 _ -)
        if not isinstance(checkpoint_id, str) or re.search(
                r"[^a-zA-Z0-9_\-]", checkpoint_id) is not None:
            raise ValueError("Invalid checkpoint id.")

        if not isinstance(checkpoint_path, (str, pathlib.Path)):
            raise ValueError("Invalid checkpoint path.")

        self.checkpoint_objects = []
        self.checkpoint_signals = []
        frm = inspect.stack()[1]
        self.calling_module = inspect.getmodule(frm[0])

        checkpoint_path = pathlib.Path(checkpoint_path) / checkpoint_id
        self.root = checkpoint_path.resolve()
        self.root.mkdir(exist_ok=True)
        self.path_signals = self.root / "signals"

        # update checkpoint counter
        self.counter = 0
        while (self.root / f"{self.counter}.checkpoint").is_file():
            self.counter += 1

        # init signals
        for signum in self.read_signals():
            self.register_signal(signum)

    def _getattr_submodule(self, obj, name, default):
        """
        Generalization of ``getattr()``.
        ``_getattr_submodule(object, "name1.sub1.sub2", None)`` will return
        attribute ``sub2`` if available otherwise ``None``.

        """
        names = name.split('.')

        for i in range(len(names) - 1):
            obj = getattr(obj, names[i], default)

        return getattr(obj, names[-1], default)

    def _setattr_submodule(self, obj, name, value):
        """
        Generalization of ``setattr()``.
        ``_setattr_submodule(object, "name1.sub1.sub2", value)`` will set
        attribute ``sub2`` to ``value``. Will raise exception if parent
        modules do not exist.

        """
        names = name.split('.')
        for i in range(len(names) - 1):
            obj = getattr(obj, names[i], None)

        if obj is None:
            raise Exception(
                f"Cannot set attribute of non existing submodules: {name}\n"
                f"Check the order you registered objects for checkpointing.")
        setattr(obj, names[-1], value)

    def __hasattr_submodule(self, obj, name):
        """
        Generalization of ``hasattr()``.
        ``__hasattr_submodule(object, "name1.sub1.sub2")`` will return
        ``True`` if submodule ``sub1`` has the attribute ``sub2``.

        """
        names = name.split('.')
        for i in range(len(names) - 1):
            obj = getattr(obj, names[i], None)

        return hasattr(obj, names[-1])

    def register(self, *args):
        """Register python objects for checkpointing.

        Parameters
        ----------
        args : list of :obj:`str`
            Names of python objects to be registered for checkpointing.

        """
        for varname in args:
            if not isinstance(varname, str):
                raise ValueError(
                    "The object that should be checkpointed is identified with its name given as a string.")

            # if not a in dir(self.calling_module):
            if not self.__hasattr_submodule(self.calling_module, varname):
                raise KeyError(
                    f"The given object '{varname}' was not found in the current scope.")

            if varname in self.checkpoint_objects:
                raise KeyError(
                    f"The given object '{varname}' is already registered for checkpointing.")

            obj = self._getattr_submodule(self.calling_module, varname, None)
            if isinstance(
                    obj, script_interface.ScriptInterfaceHelper) and not obj._so_checkpointable:
                raise TypeError(
                    f"Objects of type {type(obj)} cannot be checkpointed.")

            self.checkpoint_objects.append(varname)

    def unregister(self, *args):
        """Unregister python objects for checkpointing.

        Parameters
        ----------
        args : list of :obj:`str`
            Names of python objects to be unregistered for checkpointing.

        """
        for varname in args:
            if not isinstance(varname, str):
                raise ValueError(
                    "The object that should be checkpointed is identified with its name given as a string.")
            if varname not in self.checkpoint_objects:
                raise KeyError(
                    f"The given object '{varname}' was not registered for checkpointing yet.")

            self.checkpoint_objects.remove(varname)

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
        :obj:`bool`
            ``True`` if any checkpoints exist that match ``checkpoint_id`` and
            ``checkpoint_path`` otherwise ``False``.

        """
        return self.counter > 0

    def get_last_checkpoint_index(self):
        """
        Returns the last index of the given checkpoint id. Will raise exception
        if no checkpoints are found.

        """
        if not self.has_checkpoints():
            raise Exception(
                "No checkpoints found. Cannot return index for last checkpoint.")
        return self.counter - 1

    def save(self, checkpoint_index=None):
        """
        Saves all registered python objects in the given checkpoint directory
        using pickle.

        """
        # get attributes of registered objects
        checkpoint_data = collections.OrderedDict()
        for obj_name in self.checkpoint_objects:
            checkpoint_data[obj_name] = self._getattr_submodule(
                self.calling_module, obj_name, None)

        if checkpoint_index is None:
            checkpoint_index = self.counter

        checkpoint_file = self.root / f"{checkpoint_index}.checkpoint"
        checkpoint_file_tmp = checkpoint_file.with_suffix(".checkpoint.tmp")
        with checkpoint_file_tmp.open("wb") as f:
            pickle.dump(checkpoint_data, f, -1)
        checkpoint_file_tmp.rename(checkpoint_file)

    def load(self, checkpoint_index=None):
        """
        Loads the python objects using (c)Pickle and sets them in the calling
        module.

        Parameters
        ----------
        checkpoint_index : :obj:`int`, optional
            If not given, the last ``checkpoint_index`` will be used.

        """
        if checkpoint_index is None:
            checkpoint_index = self.get_last_checkpoint_index()

        checkpoint_file = self.root / f"{checkpoint_index}.checkpoint"
        with checkpoint_file.open("rb") as f:
            checkpoint_data = pickle.load(f)

        for key in checkpoint_data:
            self._setattr_submodule(
                self.calling_module, key, checkpoint_data[key])
            self.checkpoint_objects.append(key)

    def _signal_handler(self, signum, frame):  # pylint: disable=unused-argument
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
        if self.path_signals.is_file():
            return [int(i) for i in self.path_signals.read_text().split()]
        return []

    def _write_signal(self, signum=None):
        """Writes the given signal integer signum to the signal file.

        """
        if not utils.is_valid_type(signum, int):
            raise ValueError("Signal must be an integer number.")

        signals = self.read_signals()

        if signum not in signals:
            signals.append(signum)
            self.path_signals.write_text(" ".join(str(i) for i in signals))

    def register_signal(self, signum=None):
        """Register a signal that will trigger the signal handler.

        Parameters
        ----------
        signum : :obj:`int`
            Signal to be registered.

        """
        if not utils.is_valid_type(signum, int):
            raise ValueError("Signal must be an integer number.")

        if signum in self.checkpoint_signals:
            raise KeyError(
                f"The signal {signum} is already registered for checkpointing.")

        signal.signal(int(signum), self._signal_handler)
        self.checkpoint_signals.append(signum)
        self._write_signal(signum)
