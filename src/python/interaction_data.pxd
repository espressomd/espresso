# Handling of interactions

from espresso cimport *
cimport numpy as np
from utils cimport *

cdef extern from "interaction_data.hpp":
  ctypedef struct IA_parameters:
    double LJ_eps
    double LJ_sig
    double LJ_cut
    double LJ_shift
    double LJ_offset
    double LJ_capradius
    double LJ_min

  cdef IA_parameters *get_ia_param(int i, int j)

cdef extern from "lj.hpp":
  cdef int lennard_jones_set_params(int part_type_a, int part_type_b,
                                        double eps, double sig, double cut,
                                        double shift, double offset,
                                        double cap_radius, double min)
  cdef int ljforcecap_set_params(double ljforcecap)

ctypedef union bond_params_union:
    #* Parameters for FENE bond Potential.
    #	k - spring constant.
    #	drmax - maximal bond streching.
    #   r0 - equilibrium bond length.
    #   drmax2 - square of drmax (internal parameter). 
    
    struct fene: 
      double k
      double drmax
      double r0
      double drmax2
      double drmax2i

    #* Parameters for hyperelastic stretching_force 
    struct stretching_force:
      double r0
      double ks
    
    #* Parameters for linear stretching_force 
    struct Cstretchlin_force:
          double r0
      double kslin
    
    #* Parameters for area_force_local 
    struct area_force_local:
      double A0_l
      double ka_l
    
    #* Parameters for area_force_global 
    struct area_force_global:
      double A0_g
      double ka_g
    
    #* Parameters for bending_force 
    struct bending_force:
          double phi0
      double kb
    
    #* Parameters for volume_force 
    struct volume_force:
      double V0
      double kv
    
    struct harmonic:
      double k
      double r
      double r_cut

    #* Parameters for three body angular potential (bond-angle potentials). 
        ATTENTION: Note that there are different implementations of the bond angle
        potential which you may chose with a compiler flag in the file \ref config.hpp !
        bend - bending constant.
        phi0 - equilibrium angle (default is 180 degrees / Pi) 

    struct angle:
      double bend
      double phi0
      double cos_phi0
      double sin_phi0

    #* Parameters for three body angular potential (bond_angle_harmonic). 
    #   bend - bending constant.
    #   phi0 - equilibrium angle (default is 180 degrees / Pi) 
    struct angle_harmonic: 
      double bend
      double phi0

    #* Parameters for three body angular potential (bond_angle_cosine). 
    #   bend - bending constant.
    #   phi0 - equilibrium angle (default is 180 degrees / Pi) 
    struct angle_cosine:
      double bend
      double phi0
      double cos_phi0
      double sin_phi0

    #* Parameters for three body angular potential (bond_angle_cossquare). 
    #   bend - bending constant.
    #   phi0 - equilibrium angle (default is 180 degrees / Pi) 
    struct angle_cossquare:
      double bend
      double phi0
      double cos_phi0

   #* Parameters for four body angular potential (dihedral-angle potentials). 
    struct dihedral:
      double mult
      double bend
      double phase
    
    #* Parameters for n-body tabulated potential (n=2,3,4). 
    struct tab:
      char   *filename
      int    type
      int    npoints
      double minval
      double maxval
      double invstepsize
      double *f
      double *e
    
    #* Parameters for n-body overlapped potential (n=2,3,4). 
    struct overlap:
      char   *filename
      int    type
      double maxval
      int    noverlaps
      double *para_a
      double *para_b
      double *para_c
    
    #* Dummy parameters for -LJ Potential 
    struct subst_lj:
      double k
      double r
      double r2
    
    #*Parameters for the rigid_bond/SHAKE/RATTLE ALGORITHM
    struct rigid_bond:
      #*Length of rigid bond/Constrained Bond
      //double d
      #*Square of the length of Constrained Bond
      double d2
      #*Positional Tolerance/Accuracy value for termination of RATTLE/SHAKE iterations during position corrections
      double p_tol
      #*Velocity Tolerance/Accuracy for termination of RATTLE/SHAKE iterations during velocity corrections 
      double v_tol

    #* Parameters for three body angular potential (bond-angle potentials) that 
    #   depends on distance to wall constraint.
    #    ATTENTION: Note that there are different implementations of the bond angle
    #   potential which you may chose with a compiler flag in the file \ref config.hpp !
    #   bend - bending constant.
    #   phi0 - equilibrium angle (default is 180 degrees / Pi)
    #    dist0 - equilibrium distance (no default) 
    struct angledist:
      double bend
      double phimin
      double distmin
      double phimax
      double distmax
      double cos_phi0
      double sin_phi0

    #* Parameters for chainend angular potential with wall  
    struct endangledist:
      double bend
      double phi0
      double distmin
      double distmax



cdef extern from "interaction_data.hpp":
 ctypedef struct Bonded_ia_parameters: 
  int type
  int num
  #* union to store the different bonded interaction parameters. */
  bond_parameter_union p

 Bonded_ia_parameters* bonded_ia_params
 cdef int n_bonded_ia

cdef extern from "fene.hpp":
  int fene_set_params(int bond_type, double k, double drmax, double r0)
