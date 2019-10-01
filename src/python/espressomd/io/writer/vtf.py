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
    This fills the gap for particle ID's as required by VMD

    Parameters
    ----------
    system: espressomd.System() object
    types : :obj:`str`
            Specifies the particle types. The id mapping depends on which
            particles are going to be printed. This should be the same as
            the one used in writevsf() and writevsf().
    Returns
    -------
    dict:   A dictionary where the values are the VTF indices and the keys are the ESPresSo particle ``id``
    """

    if not hasattr(types, '__iter__'):
        types = [types]
    if types == "all":
        types = [types]
    id_to_write = []
    for p in system.part:
        for t in types:
            if p.type == t or t == "all":
                id_to_write.append(p.id)
    return dict(zip(id_to_write, range(len(id_to_write))))


def writevsf(system, fp, types='all'):
    """
    writes a VST (VTF Structure Format) to a file.
    This can be used to write the header of a VTF file.

    Parameters
    ----------
    system: espressomd.System() object
    types : :obj:`str`
            Specifies the particle types. The string 'all' will write all particles
    fp : file
               File pointer to write to.

    """

    vtf_index = vtf_pid_map(system, types)
    fp.write("unitcell {} {} {}\n".format(*(system.box_l)))

    for pid, vtf_id, in vtf_index.items():
        fp.write("atom {} radius 1 name {} type {} \n".format(vtf_id,
                                                              system.part[
                                                                  pid].type,
                                                              system.part[pid].type))
    for pid, vtf_id, in vtf_index.items():
        for b in system.part[pid].bonds:
            if (system.part[b[1]].id in vtf_index):
                fp.write("bond {}:{}\n".format(
                    vtf_id, vtf_index[system.part[b[1]].id]))


def writevcf(system, fp, types='all'):
    """
    writes a VCF (VTF Coordinate Format) to a file.
    This can be used to write a timestep to a VTF file.

    Parameters
    ----------
    system: espressomd.System() object
    types : :obj:`str`
            Specifies the particle types. The string 'all' will write all particles
    fp : file
               File pointer to write to.

    """
    vtf_index = vtf_pid_map(system, types)
    fp.write("\ntimestep indexed\n")
    for pid, vtf_id, in vtf_index.items():
        fp.write("{} {} {} {}\n".format(vtf_id, *(system.part[pid].pos)))
