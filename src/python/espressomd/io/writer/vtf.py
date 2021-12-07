# Copyright (C) 2010-2019 The ESPResSo project
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


def vtf_pid_map(system, types='all'):
    """
    Generates a VTF particle index map to ESPResSo ``id``.
    This fills the gap for particle ID's as required by VMD.

    Parameters
    ----------
    system: :obj:`espressomd.system.System`
    types : :obj:`str`
        Specifies the particle types. The id mapping depends on which
        particles are going to be printed. This should be the same as
        the one used in :func:`writevsf()` and :func:`writevcf()`.
    Returns
    -------
    dict:
        A dictionary where the values are the VTF indices and the keys
        are the ESPresSo particle ``id``.
    """

    if not hasattr(types, '__iter__'):
        types = [types]
    if types == "all":
        types = [types]
    id_to_write = []
    for p in system.part:
        for t in types:
            if t in (p.type, "all"):
                id_to_write.append(p.id)
    return dict(zip(id_to_write, range(len(id_to_write))))


def writevsf(system, fp, types='all'):
    """
    writes a VST (VTF Structure Format) to a file.
    This can be used to write the header of a VTF file.

    Parameters
    ----------
    system: :obj:`espressomd.system.System`
    types : :obj:`str`
        Specifies the particle types. The string 'all' will write all particles
    fp : file
        File pointer to write to.

    """

    vtf_index = vtf_pid_map(system, types)
    fp.write(f"unitcell {' '.join(map(str, system.box_l))}\n")

    for pid, vtf_id, in vtf_index.items():
        partcl = system.part.by_id(pid)
        fp.write(
            f"atom {vtf_id} radius 1 name {partcl.type} type {partcl.type} \n")
    for pid, vtf_id, in vtf_index.items():
        for b in system.part.by_id(pid).bonds:
            if system.part.by_id(b[1]).id in vtf_index:
                fp.write(
                    f"bond {vtf_id}:{vtf_index[system.part.by_id(b[1]).id]}\n")


def writevcf(system, fp, types='all'):
    """
    writes a VCF (VTF Coordinate Format) to a file.
    This can be used to write a timestep to a VTF file.

    Parameters
    ----------
    system: :obj:`espressomd.system.System`
    types : :obj:`str`
        Specifies the particle types. The string 'all' will write all particles
    fp : file
        File pointer to write to.

    """
    vtf_index = vtf_pid_map(system, types)
    fp.write("\ntimestep indexed\n")
    for pid, vtf_id, in vtf_index.items():
        fp.write(f"{vtf_id} {' '.join(map(str, system.part.by_id(pid).pos))}\n")
