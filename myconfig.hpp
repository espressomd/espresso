/*
    This is a sample for the file @myconfig@.

    Create a copy of this file in any of these locations: 
        @myconfig@ 
	/home/itvgi/PHD/CODE/TPB_NEW/three-particle-binding_new/@myconfig@

    Then uncomment any of the following lines to activate the
    corresponding feature of ESPResSo. It is recommended to activate
    only those features that you actually need to optimize the
    performance of ESPResSo for your problem.

    To access the information on the compilation status of the code
    you are working with in your ESPResSo Tcl-script, use the
    corresponding \ref tcl_features "Tcl-commands".
*/

/* Generic features */
#define PARTIAL_PERIODIC
//#define EXTERNAL_FORCES
//#define CONSTRAINTS
//#define MASS
//#define EXCLUSIONS
//#define COMFORCE
//#define COMFIXED
//#define MOLFORCES
//#define BOND_CONSTRAINT
//#define MODES
//#define BOND_VIRTUAL
//#define LANGEVIN_PER_PARTICLE
#define COLLISION_DETECTION
//#define METADYNAMICS
//#define NEMD
//#define NPT
//#define GHMC
//#define CATALYTIC_REACTIONS

//#define HARMONICFORCE

/* Rotation */
//#define ROTATION
//#define ROTATIONAL_INERTIA
//#define ROTATION_PER_PARTICLE

/* Adress */
//#define MOL_CUT

/* Electrostatics */
//#define ELECTROSTATICS
//#define INTER_RF

/* Magnetostatics */
//#define DIPOLES

/* Virtual sites features */
//#define VIRTUAL_SITES_COM
#define VIRTUAL_SITES_RELATIVE
//#define VIRTUAL_SITES_NO_VELOCITY
//#define VIRTUAL_SITES_THERMOSTAT
//#define THERMOSTAT_IGNORE_NON_VIRTUAL

/* DPD features */
//#define DPD
//#define TRANS_DPD
//#define DPD_MASS_RED
//#define DPD_MASS_LIN
//#define TUNABLE_SLIP

/* Alternative implementation of DPD */
//#define INTER_DPD

/* Lattice-Boltzmann features */
//#define LB
//#define LB_GPU
//#define LB_BOUNDARIES
//#define LB_BOUNDARIES_GPU
//#define LB_ELECTROHYDRODYNAMICS
//#define ELECTROKINETICS
//#define EK_BOUNDARIES
//#define EK_REACTION
//#define SHANCHEN

/* Interaction features */
//#define TABULATED
#define LENNARD_JONES
//#define LJ_WARN_WHEN_CLOSE
//#define LENNARD_JONES_GENERIC
//#define LJCOS
//#define LJCOS2
//#define LJ_ANGLE
//#define LJGEN_SOFTCORE
//#define GAY_BERNE
//#define SMOOTH_STEP
//#define HERTZIAN
//#define GAUSSIAN
//#define BMHTF_NACL
//#define MORSE
//#define BUCKINGHAM
//#define SOFT_SPHERE
//#define HAT
//#define OVERLAPPED

/* Fluid-Structure Interactions (object in fluid) */
//#define AREA_FORCE_GLOBAL
//#define VOLUME_FORCE

#define BOND_ANGLE
#define COLLISION_USE_BROKEN_PARALLELIZATION

//#define BOND_ANGLEDIST
//#define BOND_ANGLEDIST_HARMONIC

//#define BOND_ENDANGLEDIST
//#define BOND_ENDANGLEDIST_HARMONIC

/* Obsolete features. To be removed. */
/* Old Bond angle */
/* Note: Activate ONLY ONE bonded angle potential out of the following! */
//#define BOND_ANGLE_HARMONIC
//#define BOND_ANGLE_COSINE
//#define BOND_ANGLE_COSSQUARE

/* Strange features. Use only if you know what you are doing! */
/* activate the old dihedral form */
//#define OLD_DIHEDRAL
/* turn off nonbonded interactions within molecules */
//#define NO_INTRA_NB


/* Debugging */
//#define ADDITIONAL_CHECKS
//#define ASYNC_BARRIER

//#define COMM_DEBUG
//#define EVENT_DEBUG
//#define INTEG_DEBUG
//#define CELL_DEBUG
//#define GHOST_DEBUG
//#define LATTICE_DEBUG
//#define HALO_DEBUG
//#define GRID_DEBUG
//#define VERLET_DEBUG
//#define PARTICLE_DEBUG
//#define P3M_DEBUG
//#define FFT_DEBUG
//#define RANDOM_DEBUG
//#define FORCE_DEBUG
//#define THERMO_DEBUG
//#define LJ_DEBUG
//#define MORSE_DEBUG
//#define ESR_DEBUG
//#define ESK_DEBUG
//#define FENE_DEBUG
//#define GHOST_FORCE_DEBUG
//#define STAT_DEBUG
//#define POLY_DEBUG
//#define MOLFORCES_DEBUG
//#define PTENSOR_DEBUG
//#define MEM_DEBUG
//#define MAGGS_DEBUG
//#define LB_DEBUG
//#define VIRTUAL_SITES_DEBUG
//#define GHOST_FLAG 
//#define MPI_CORE
//#define FORCE_CORE

/* Single particle debugging */
//#define ONEPART_DEBUG
// which particle id to debug
//#define ONEPART_DEBUG_ID 13

