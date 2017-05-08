def writevsf(system, fp, types='all'):
    """
    writes a VST (VTF Structure Format) to a file.
    This can be used to write the header of a VTF file.

    Parameters
    ----------
    system: espressomd.System() object
    types : str
            Specifies the particle types. The string 'all' will write all particles
    fp : file
               File pointer to write to.

    """

    if not hasattr(types, '__iter__'):
        types = [types]
    if types == "all":
        types = [types]
        
    fp.write("unitcell {} {} {}\n".format(*(system.box_l)))
    for p in system.part:
        for t in types:
            if (p.type == t or t == "all"):
                fp.write("atom {} radius 1 name {} type {} \n".format(p.id, p.type, p.type))

    for p in system.part:
        for t in types:
            if (p.type == t or t == "all"):
                for b in p.bonds:
                    if (system.part[b[1]].type == t or t == "all"):
                        fp.write("bond {}:{}\n".format(p.id,system.part[b[1]].id))


def writevcf(system, fp, types='all'):
    """
    writes a VCF (VTF Coordinate Format) to a file.
    This can be used to write a stimestep to a VTF file.
    
    Parameters
    ----------
    system: espressomd.System() object
    types : str
            Specifies the particle types. The string 'all' will write all particles
    fp : file
               File pointer to write to.

    """

    if not hasattr(types, '__iter__'):
        types = [types]
    if types == "all":
        types = [types]

    fp.write("\ntimestep indexed\n")
    for p in system.part:
        for t in types:
            if (p.type == t or t == "all"):
                fp.write("{} {} {} {}\n".format(p.id, *(p.pos)))
