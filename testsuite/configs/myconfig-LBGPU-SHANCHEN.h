/* minimal config to run LB_GPU testcase */
#define LB_GPU
#define SHANCHEN
// not yet possible, if so, merge with myconfig-LBGPU
//#define LB_BOUNDARIES_GPU

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
#define MODES
#define BOND_VIRTUAL
#define COLLISION_DETECTION
#define LANGEVIN_PER_PARTICLE
#define ROTATION_PER_PARTICLE
#define CATALYTIC_REACTIONS
#define GRANDCANONCIAL

#define NEMD
#define NPT 
#define GHMC

#define AREA_FORCE_GLOBAL   
#define VOLUME_FORCE   

#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
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
