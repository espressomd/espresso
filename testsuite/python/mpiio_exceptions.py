#
# Copyright (C) 2022 The ESPResSo project
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

import espressomd
import espressomd.io

import unittest as ut
import os
import tempfile


class MPIIOMockGenerator:
    """Mock MPI-IO files."""

    def __init__(self, root):
        self.root = root
        self.counter = 0

    def _create_file(self, fn, mode, content=b''):
        with open(os.open(fn, os.O_CREAT | os.O_WRONLY, mode), 'wb') as f:
            f.write(content)

    def create(self, *suffixes, read_only=True, from_ref=None):
        mode = 0o444 if read_only else 0o777
        filepath = os.path.join(self.root, f'testdata.mpiio.{self.counter}')
        filepath_derived = []
        for suffix in suffixes:
            derived = f'{filepath}.{suffix}'
            content = b''
            if from_ref is not None:
                with open(f'{from_ref}.{suffix}', 'rb') as f:
                    content = f.read()
            self._create_file(derived, mode, content)
            filepath_derived.append(derived)
        if len(filepath_derived) == 1:
            filepath_derived = filepath_derived[0]
        self.counter += 1
        return filepath, filepath_derived


class MPIIOTest(ut.TestCase):

    """
    Test class for the MPI-IO core functionality.
    Check for exceptions when data cannot be read or written.
    With 1 MPI rank, fatal errors are just exceptions.
    """
    system = espressomd.system.System(box_l=[1, 1, 1])
    n_nodes = system.cell_system.get_state()["n_nodes"]

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    @ut.skipIf(n_nodes != 1, "only works on 1 MPI rank")
    def test_exceptions(self):
        generator = MPIIOMockGenerator(self.temp_dir.name)
        mpiio = espressomd.io.mpiio.Mpiio()

        # generate reference data
        self.system.part.add(pos=[0, 0, 0])
        path_ref = generator.create()[0]
        mpiio.write(path_ref, types=True)
        self.system.part.clear()

        # check reference data is valid
        mpiio.read(path_ref, types=True)
        self.system.part.clear()

        # exception when the metadata cannot be written
        path, fn = generator.create('head', read_only=True)
        with self.assertRaisesRegex(RuntimeError, f'Could not open file "{fn}"'):
            mpiio.write(path, types=True)

        # exception when the payload cannot be written
        path, fn = generator.create('pref', read_only=True)
        with self.assertRaisesRegex(RuntimeError, f'Could not open file "{fn}"'):
            mpiio.write(path, types=True)

        # exception when calculating the size of a non-existent file
        path, _ = generator.create(read_only=True)
        fn = f'{path}.pref'
        with self.assertRaisesRegex(RuntimeError, f'Could not get file size of "{fn}"'):
            mpiio.read(path, types=True)

        # exception when the MPI world size differs for reading and writing
        # (empty .pref file -> data was written with MPI world size of 0)
        path, _ = generator.create('id', 'pref', read_only=True)
        with self.assertRaisesRegex(RuntimeError, f'Trying to read a file with a different COMM size than at point of writing'):
            mpiio.read(path, types=True)

        # exception when the particle types don't exist
        path, _ = generator.create(
            'pref', 'id', 'head', read_only=False, from_ref=path_ref)
        fn = f'{path}.type'
        with self.assertRaisesRegex(RuntimeError, f'Could not open file "{fn}"'):
            mpiio.read(path, types=True)

        # exception when the metadata doesn't exist
        path, _ = generator.create(
            'id', 'pref', read_only=False, from_ref=path_ref)
        fn = f'{path}.head'
        with self.assertRaisesRegex(RuntimeError, f'Could not open file "{fn}"'):
            mpiio.read(path, types=True)

        # exception when the metadata is empty
        with open(fn, 'wb'):
            pass
        with self.assertRaisesRegex(RuntimeError, f'Could not read file "{fn}"'):
            mpiio.read(path, types=True)

        # exception when reading data that was not written to disk
        # (empty .pref file -> data was written with MPI world size of 0)
        path, _ = generator.create(
            'id', 'pref', 'head', read_only=False, from_ref=path_ref)
        with self.assertRaisesRegex(RuntimeError, f'Requesting to read fields which were not dumped'):
            mpiio.read(path, types=True, bonds=True)


if __name__ == '__main__':
    ut.main()
