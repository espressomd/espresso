/* maximal set of features usable at the same time plus all debug switches */
/* Do not run the testsuite with this set, only compile it. */
#define PARTIAL_PERIODIC
#define ELECTROSTATICS
#define DIPOLES
#define ROTATION
#define ROTATIONAL_INERTIA
#define PARTICLE_ANISOTROPY
#define EXTERNAL_FORCES

#define MASS
#define EXCLUSIONS
#define MOLFORCES


#define COLLISION_DETECTION
#define LANGEVIN_PER_PARTICLE
#define SWIMMER_REACTIONS

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
#define MMM1D_GPU
#endif

#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define LJGEN_SOFTCORE
#define LJCOS
#define LJCOS2
#define GAUSSIAN
#define HAT
#define GAY_BERNE
#define SMOOTH_STEP
#define HERTZIAN
#define BMHTF_NACL
#define MORSE
#define BUCKINGHAM
#define SOFT_SPHERE
#define INTER_RF
#define OVERLAPPED
#define THOLE

#define BOND_ANGLE
#define BOND_ANGLEDIST
#define BOND_ANGLEDIST_HARMONIC

#define VIRTUAL_SITES_RELATIVE

#define EXPERIMENTAL_FEATURES
// DEBUG Switches
#define LJ_WARN_WHEN_CLOSE

#define COMM_DEBUG
#define EVENT_DEBUG
#define INTEG_DEBUG
#define CELL_DEBUG
#define GHOST_DEBUG
#define LATTICE_DEBUG
#define HALO_DEBUG
#define GRID_DEBUG
#define VERLET_DEBUG
#define PARTICLE_DEBUG
#define P3M_DEBUG
#define FFT_DEBUG
#define RANDOM_DEBUG
#define FORCE_DEBUG
#define THERMO_DEBUG
#define LJ_DEBUG
#define MORSE_DEBUG
#define ESR_DEBUG
#define ESK_DEBUG
#define FENE_DEBUG
#define GHOST_FORCE_DEBUG
#define STAT_DEBUG
#define POLY_DEBUG
#define MOLFORCES_DEBUG
#define PTENSOR_DEBUG
#define MAGGS_DEBUG
#define LB_DEBUG
#define VIRTUAL_SITES_DEBUG
#define LE_DEBUG
#ifdef CUDA
#define CUDA_DEBUG
#endif
#define ESIF_DEBUG

#define ONEPART_DEBUG
