/* test configuration without relative/com virtual sites for immersed boundaries */
#define IMMERSED_BOUNDARY
#define PARTIAL_PERIODIC
#define ELECTROSTATICS
#define DIPOLES
#define ROTATION
#define ROTATIONAL_INERTIA
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

#define BOND_ANGLE
#define BOND_ANGLEDIST
#define BOND_ANGLEDIST_HARMONIC
