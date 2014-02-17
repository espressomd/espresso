/* maximal set of features usable at the same time */
#define PARTIAL_PERIODIC
#define ELECTROSTATICS
#define DIPOLES
#define ROTATION
#define ROTATIONAL_INERTIA
#define EXTERNAL_FORCES
#define CONSTRAINTS
#define MASS
#define EXCLUSIONS
#define COMFORCE
#define COMFIXED
#define MOLFORCES

#ifdef FFTW
#define MODES
#endif

#define BOND_VIRTUAL
#define COLLISION_DETECTION
#define LANGEVIN_PER_PARTICLE
#define ROTATION_PER_PARTICLE
#define CATALYTIC_REACTIONS

#define NEMD
#define NPT 
#define GHMC

#define LB
#define LB_BOUNDARIES
#define LB_ELECTROHYDRODYNAMICS

#ifdef CUDA
#define LB_GPU
#define LB_BOUNDARIES_GPU
#define ELECTROKINETICS
#define EK_BOUNDARIES
#define EK_REACTION
#endif

#define AREA_FORCE_GLOBAL   
#define VOLUME_FORCE   

#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define LJGEN_SOFTCORE
#define LJCOS
#define LJCOS2
#define GAUSSIAN
#define HAT
#define LJ_ANGLE
#define GAY_BERNE
#define SMOOTH_STEP
#define HERTZIAN
#define BMHTF_NACL
#define MORSE
#define BUCKINGHAM
#define SOFT_SPHERE
#define INTER_RF
#define OVERLAPPED

#define BOND_ANGLE
#define BOND_ANGLEDIST
#define BOND_ANGLEDIST_HARMONIC
#define BOND_ENDANGLEDIST
#define BOND_ENDANGLEDIST_HARMONIC

#define VIRTUAL_SITES_RELATIVE
