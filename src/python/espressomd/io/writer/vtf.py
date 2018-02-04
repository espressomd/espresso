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
    dict:   A dictionary where the values are the VTF indicies and the keys are the ESPresSo particle ``id``
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

    for pid, vtf_id, in vtf_index.iteritems():
        fp.write("atom {} radius 1 name {} type {} \n".format(vtf_id,
                                                              system.part[pid].type,
                                                              system.part[pid].type))
    for pid, vtf_id, in vtf_index.iteritems():
        for b in system.part[pid].bonds:
            if (system.part[b[1]].id in vtf_index ):
                fp.write("bond {}:{}\n".format(vtf_id, vtf_index[system.part[b[1]].id]))


def writevcf(system, fp, types='all'):
    """
    writes a VCF (VTF Coordinate Format) to a file.
    This can be used to write a stimestep to a VTF file.
    
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
    for pid, vtf_id, in vtf_index.iteritems():
        fp.write("{} {} {} {}\n".format(vtf_id, *(system.part[pid].pos)))
