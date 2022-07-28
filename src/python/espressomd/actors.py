# Copyright (C) 2010-2022 The ESPResSo project
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

from . import highlander
from . import script_interface


class Actors:

    """
    Container for actor objects.
    """

    active_actors = []

    def __del__(self):
        self.clear()

    def __getstate__(self):
        return self.active_actors

    def __setstate__(self, active_actors):
        self.active_actors[:] = []
        for actor in active_actors:
            self.active_actors.append(actor)
            actor._activate()

    def add(self, actor):
        """
        Parameters
        ----------
        actor :
            Actor to add to this container.

        """
        if actor in Actors.active_actors:
            raise highlander.ThereCanOnlyBeOne(actor)

        if isinstance(actor, script_interface.ScriptInterfaceHelper):
            actor._activate()

        self.active_actors.append(actor)

        if not isinstance(actor, script_interface.ScriptInterfaceHelper):
            actor._activate()

    def remove(self, actor):
        """
        Parameters
        ----------
        actor :
            Actor to remove from this container.

        """
        if actor not in self.active_actors:
            raise Exception("Actor is not active")
        actor._deactivate()
        self.active_actors.remove(actor)

    def clear(self):
        """Remove all actors."""
        # The order in which actors are removed matters. Some actors set up
        # global bitfields that activate sanity checks on the MD cellsystem,
        # and reset these bitfields when removed. Actors need to be removed
        # in the reverse order they were inserted to guarantee pre-conditions
        # and post-conditions are always met.
        while len(self.active_actors):
            self.remove(self.active_actors[-1])

    def __str__(self):
        return str(self.active_actors)

    def __getitem__(self, key):
        return self.active_actors[key]

    def __len__(self):
        return len(self.active_actors)

    def __iter__(self):
        for actor in self.active_actors:
            yield actor

    def __delitem__(self, idx):
        actor = self[idx]
        self.remove(actor)
