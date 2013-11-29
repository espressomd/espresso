import sys

FEATURES=[ \
"PARTIAL_PERIODIC", \
"MAGNETOSTATICS", \
"ELECTROSTATICS", \
"ROTATION", \
"ROTATIONAL_INERTIA", \
"DIPOLES", \
"EXTERNAL_FORCES", \
"CONSTRAINTS", \
"MASS", \
"EXCLUSIONS", \
"COMFORCE", \
"COMFIXED", \
"MOLFORCES", \
"BOND_CONSTRAINT", \
"MODES", \
"BOND_VIRTUAL", \
"LANGEVIN_PER_PARTICLE", \
"ADRESS", \
"METADYNAMICS", \
"OVERLAPPED", \
"VIRTUAL_SITES_COM", \
"VIRTUAL_SITES_RELATIVE", \
"VIRTUAL_SITES_NO_VELOCITY", \
"VIRTUAL_SITES_THERMOSTAT", \
"THERMOSTAT_IGNORE_NON_VIRTUAL", \
"NEMD", \
"NPT ", \
"DPD", \
"TRANS_DPD", \
"INTER_DPD", \
"DPD_MASS_RED", \
"DPD_MASS_LIN", \
"LB", \
"LB_GPU", \
"LB_BOUNDARIES", \
"LB_BOUNDARIES_GPU", \
"LB_ELECTROHYDRODYNAMICS", \
"TABULATED", \
"LENNARD_JONES", \
"LJ_WARN_WHEN_CLOSE", \
"LENNARD_JONES_GENERIC", \
"LJCOS", \
"LJCOS2", \
"LJ_ANGLE", \
"GAY_BERNE", \
"SMOOTH_STEP", \
"HERTZIAN", \
"BMHTF_NACL", \
"MORSE", \
"BUCKINGHAM", \
"SOFT_SPHERE", \
"INTER_RF", \
"MOL_CUT", \
"TUNABLE_SLIP", \
"NO_INTRA_NB", \
"BOND_ANGLE_HARMONIC", \
"BOND_ANGLE_COSINE", \
"BOND_ANGLE_COSSQUARE", \
"BOND_ANGLEDIST", \
"BOND_ENDANGLEDIST", \
"OLD_DIHEDRAL", \
"P3M_BRILLOUIN ", \
"P3M_MAX_MESH", \
"USE_ERFC_APPROXIMATION ", \
"ROUND_ERROR_PREC", \
"TINY_SIN_VALUE", \
"TINY_COS_VALUE", \
"TINY_LENGTH_VALUE", \
"SHAKE_MAX_ITERATIONS", \
"ADDITIONAL_CHECKS", \
"COMM_DEBUG", \
"EVENT_DEBUG", \
"INTEG_DEBUG", \
"CELL_DEBUG", \
"GHOST_DEBUG", \
"LATTICE_DEBUG", \
"HALO_DEBUG", \
"GRID_DEBUG", \
"VERLET_DEBUG", \
"PARTICLE_DEBUG", \
"P3M_DEBUG", \
"EWALD_DEBUG", \
"FFT_DEBUG", \
"RANDOM_DEBUG", \
"FORCE_DEBUG", \
"THERMO_DEBUG ", \
"LJ_DEBUG", \
"MORSE_DEBUG", \
"ESR_DEBUG", \
"ESK_DEBUG", \
"FENE_DEBUG", \
"GHOST_FORCE_DEBUG", \
"ONEPART_DEBUG", \
"STAT_DEBUG", \
"POLY_DEBUG", \
"MOLFORCES_DEBUG", \
"PTENSOR_DEBUG", \
"MEM_DEBUG", \
"MAGGS_DEBUG", \
"LB_DEBUG", \
"VIRTUAL_SITES_DEBUG", \
"ASYNC_BARRIER", \
"MPI_CORE", \
"FORCE_CORE", \
"OLD_RW_VERSION" \
]
import re
def feature_value(myconfig, f):
    without_c_comments=re.sub("/\*.*?\*/", " ", myconfig, flags=re.DOTALL)
    without_cpp_comments=re.sub("//.*", " ", without_c_comments)
    # The following regex requires the string to have at least a space or line break after the constant name
    m=re.match(".*(#define)\s*("+f+")\s+(\.*)", without_cpp_comments + "\n", re.DOTALL)
    print f,m
    if m:
        if m.group(3) != "":
            return m.group(3)
        else:
            return "1"
    else:
        return "0"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("Usage: make_myconfig.py <myconfig.h>")
    else:
        myconfigfilename = sys.argv[1]
    myconfig = open(myconfigfilename).read()
    pxifilename="myconfig.pxi"
    pxifile=open(pxifilename, "w")
    for f in FEATURES:
        value = feature_value(myconfig, f)
        pxifile.write("DEF " + f + " = " + value + "\n") 

  
  
