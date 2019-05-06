# Copyright (C) 2010-2018 The ESPResSo project
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
"""
This modules allows to expose ESPREsSo's coordinates and particle attributes
to MDAnalysis without need to save information to files.

The main class is :class:`Stream`, which is used to initialize the stream of
data to MDAnalysis' readers. These are the topology reader :class:`ESPParser`
and the coordinates reader :class:`ESPReader`.

A minimal working example is the following:

>>> # imports
>>> import espressomd
>>> from espressomd import MDA_ESP
>>> import MDAnalysis as mda

>>> # system setup
>>> system = espressomd.System()
>>> system.time_step = 1.
>>> system.cell_system.skin = 1.
>>> system.box_l = [10.,10.,10.]
>>> system.part.add(id=0,pos=[1.,2.,3.])

>>> # set up the stream
>>> eos = MDA_ESP.Stream(system)
>>> # feed Universe with a topology and with coordinates
>>> u = mda.Universe(eos.topology,eos.trajectory)
>>> print u
<Universe with 1 atoms>

"""

try:
    import cStringIO as StringIO
    StringIO = StringIO.StringIO
except ImportError:
    from io import StringIO

import numpy as np
import MDAnalysis

from distutils.version import LooseVersion

from MDAnalysis.lib import util
from MDAnalysis.coordinates.core import triclinic_box
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.coordinates import base
from MDAnalysis.coordinates.base import SingleFrameReaderBase
from MDAnalysis.core.topology import Topology

from MDAnalysis.core.topologyattrs import (
    Atomnames, Atomids, Atomtypes, Masses,
    Resids, Resnums, Segids, Resnames, AltLocs,
    ICodes, Occupancies, Tempfactors, Charges
)


class Stream(object):

    """
    Create an object that provides a MDAnalysis topology and a coordinate reader

    >>> eos = MDA_ESP.Stream(system)
    >>> u = mda.Universe(eos.topology,eos.trajectory)

    Parameters
    ----------

    system : :obj:`espressomd.system.System`
    """

    def __init__(self, system):
        self.topology = ESPParser(None, espresso=system).parse()
        self.system = system

    @property
    def trajectory(self):
        """
        Particles' coordinates at the current time

        Returns
        -------

        stream : :class:`MDAnalysis.lib.util.NamedStream`
            A stream in the format that can be parsed by :class:`ESPReader`

        """
        # time
        _xyz = str(self.system.time) + '\n'
        # number of particles
        _xyz += str(len(self.system.part)) + '\n'
        # box edges
        _xyz += str(self.system.box_l) + '\n'
        # configuration
        for _p in self.system.part:
            _xyz += str(_p.pos) + '\n'
        for _p in self.system.part:
            _xyz += str(_p.v) + '\n'
        for _p in self.system.part:
            _xyz += str(_p.f) + '\n'
        return NamedStream(StringIO(_xyz), "__.ESP")


class ESPParser(TopologyReaderBase):

    """
    An MDAnalysis reader of espresso's topology

    """
    format = 'ESP'

    def __init__(self, filename, **kwargs):
        self.kwargs = kwargs

    def parse(self):
        """
        Access ESPResSo data and return the topology object

        Returns
        -------

        top : :class:`MDAnalysis.core.topology.Topology`
            a topology object

        """
        espresso = self.kwargs['espresso']

        names = []
        atomtypes = []
        masses = []
        charges = []

        for p in espresso.part:
            names.append("A" + repr(p.type))
            atomtypes.append("T" + repr(p.type))
            masses.append(p.mass)
            charges.append(p.q)
        natoms = len(espresso.part)
        attrs = [Atomnames(np.array(names, dtype=object)),
                 Atomids(np.arange(natoms) + 1),
                 Atomtypes(np.array(atomtypes, dtype=object)),
                 Masses(masses),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['System'], dtype=object)),
                 AltLocs(np.array([' '] * natoms, dtype=object)),
                 Resnames(np.array(['R'], dtype=object)),
                 Occupancies(np.zeros(natoms)),
                 Tempfactors(np.zeros(natoms)),
                 ICodes(np.array([' '], dtype=object)),
                 Charges(np.array(charges)),
                 ]

        top = Topology(natoms, 1, 1, attrs=attrs)

        return top


class Timestep(base.Timestep):
    _ts_order_x = [0, 3, 4]
    _ts_order_y = [5, 1, 6]
    _ts_order_z = [7, 8, 2]

    def _init_unitcell(self):
        return np.zeros(9, dtype=np.float32)

    @property
    def dimensions(self):
        # This information now stored as _ts_order_x/y/z to keep DRY
        x = self._unitcell[self._ts_order_x]
        y = self._unitcell[self._ts_order_y]
        z = self._unitcell[self._ts_order_z]
            # this ordering is correct! (checked it, OB)
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        x, y, z = triclinic_vectors(box)
        np.put(self._unitcell, self._ts_order_x, x)
        np.put(self._unitcell, self._ts_order_y, y)


class ESPReader(SingleFrameReaderBase):

    """
    An MDAnalysis single frame reader for the stream provided by Stream()

    """
    format = 'ESP'
    units = {'time': None, 'length': 'nm', 'velocity': 'nm/ps'}
    _Timestep = Timestep

    def _read_first_frame(self):
        with util.openany(self.filename, 'rt') as espfile:
            n_atoms = 1
            for pos, line in enumerate(espfile, start=-3):
                if (pos == -3):
                    time = float(line[1:-1])
                elif(pos == -2):
                    n_atoms = int(line)
                    self.n_atoms = n_atoms
                    positions = np.zeros(
                        self.n_atoms * 3, dtype=np.float32).reshape(self.n_atoms, 3)
                    velocities = np.zeros(
                        self.n_atoms * 3, dtype=np.float32).reshape(self.n_atoms, 3)
                    forces = np.zeros(
                        self.n_atoms * 3, dtype=np.float32).reshape(self.n_atoms, 3)
                    self.ts = ts = self._Timestep(
                        self.n_atoms, **self._ts_kwargs)
                    self.ts.time = time
                elif(pos == -1):
                    self.ts._unitcell[:3] = np.array(
                        list(map(float, line[1:-2].split())))
                elif(pos < n_atoms):
                    positions[pos] = np.array(
                        list(map(float, line[1:-2].split())))
                elif(pos < 2 * n_atoms):
                    velocities[pos - n_atoms] = np.array(
                        list(map(float, line[1:-2].split())))
                else:
                    forces[pos - 2 * n_atoms] = np.array(
                        list(map(float, line[1:-2].split())))

            ts.positions = np.copy(positions)
            ts.velocities = np.copy(velocities)
            ts.forces = np.copy(forces)
