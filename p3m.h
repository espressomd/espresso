// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef P3M_H 
#define P3M_H
/** \file p3m.h   main header-file for P3M algorithms intended to deal with long-range forces.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  Further reading: 
 *  <ul>
 *  <li> P.P. Ewald,
 *       <i>Die Berechnung optischer und elektrostatischer Gitterpotentiale</i>,
 *       Ann. Phys. (64) 253-287, 1921
 *  <li> R. W. Hockney and J. W. Eastwood, 
 *       <i>Computer Simulation Using Particles</i>,
 *       IOP, London, 1988
 *  <li> M. Deserno and C. Holm,
 *       <i>How to mesh up {E}wald sums. I. + II.</i>,
 *       J. Chem. Phys. (109) 7678, 1998; (109) 7694, 1998
 *  <li> M. Deserno, C. Holm and H. J. Limbach,
 *       <i>How to mesh up {E}wald sums. </i>,
 *       in Molecular Dynamics on Parallel Computers,
 *       Ed. R. Esser et al., World Scientific, Singapore, 2000
 *  <li> M. Deserno,
 *       <i>Counterion condensation for rigid linear polyelectrolytes</i>,
 *       PhdThesis, Universit{\"a}t Mainz, 2000
 *  <li> J.J. Cerda, P3M for dipolar interactions. J. Chem. Phys, 129, xxx ,(2008).
 *  </ul>
 *
 *  For more information about the p3m algorithm,
  *  For more information about the p3m algorithm,
 *  see \ref p3m.h "p3m.h"
 *  see \ref p3m-charges.c  "p3m-charges.c"
 *  see \ref p3m-charges.h  "p3m-charges.h"
 *  see \ref p3m-dipoles.c  "p3m-dipoles.c"
 *  see \ref p3m-dipoles.h  "p3m-dipoles.h"
 *  see \ref p3m-assignment.c  "p3m-assignment.c"
 */

#include "utils.h"
#include "interaction_data.h"
#include "integrate.h"

#ifdef ELP3M

/* defined only within this header file, for checking that system (electro/magnetostatic)
   specific files are only included from here. */
#define P3M_H_CURRENT

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
#ifdef ELECTROSTATICS
  /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double alpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double r_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    mesh[3];
  /** offset of the first mesh point (lower left 
      corner) from the coordinate origin ([0,1[). */
  double mesh_off[3];
  /** charge assignment order ([0,7]). */
  int    cao;
  /** number of interpolation points for charge assignment function */
  int    inter;
  /** Accuracy of the actual parameter set. */
  double accuracy;

  /** epsilon of the "surrounding dielectric". */
  double epsilon;
  /** Cutoff for charge assignment. */
  double cao_cut[3];
  /** mesh constant. */
  double a[3];
  /** inverse mesh constant. */
  double ai[3];
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double alpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double r_cut;
  /** full size of the interpolated assignment function */
  int inter2;
  /** number of points unto which a single charge is interpolated, i.e. p3m.cao^3 */
  int cao3;
  /** additional points around the charge assignment mesh, for method like dielectric ELC
      creating virtual charges. */
  double additional_mesh[3];
#endif  
#ifdef MAGNETOSTATICS 
    /** Ewald splitting parameter (0<alpha<1), rescaled to alpha_L = alpha * box_l. */
  double Dalpha_L;
  /** Cutoff radius for real space electrostatics (>0), rescaled to r_cut_iL = r_cut * box_l_i. */
  double Dr_cut_iL;
  /** number of mesh points per coordinate direction (>0). */
  int    Dmesh[3];
  /** offset of the first mesh point (lower left 
      corner) from the coordinate origin ([0,1[). */
  double Dmesh_off[3];
  /** charge assignment order ([0,7]). */
  int    Dcao;
  /** number of interpolation points for charge assignment function */
  int    Dinter;
  /** Accuracy of the actual parameter set. */
  double Daccuracy;

  /** epsilon of the "surrounding dielectric". */
  double Depsilon;
  /** Cutoff for charge assignment. */
  double Dcao_cut[3];
  /** mesh constant. */
  double Da[3];
  /** inverse mesh constant. */
  double Dai[3];
  /** unscaled \ref alpha_L for use with fast inline functions only */
  double Dalpha;
  /** unscaled \ref r_cut_iL for use with fast inline functions only */
  double Dr_cut;
  /** full size of the interpolated assignment function */
  int Dinter2;
  /** number of points unto which a single charge is interpolated, i.e. p3m.Dcao^3 */
  int Dcao3;
  /** additional points around the charge assignment mesh, for method like dielectric ELC
      creating virtual charges. */
  double Dadditional_mesh[3];
#endif
} p3m_struct;

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  /** dimension (size) of local mesh. */
  int dim[3];
  /** number of local mesh points. */
  int size;
  /** index of lower left corner of the 
      local mesh in the global mesh. */
  int ld_ind[3];
  /** position of the first local mesh point. */
  double ld_pos[3];
  /** dimension of mesh inside node domain. */
  int inner[3];
  /** inner left down grid point */
  int in_ld[3];
  /** inner up right grid point + (1,1,1)*/
  int in_ur[3];
  /** number of margin mesh points. */
  int margin[6];
  /** number of margin mesh points from neighbour nodes */
  int r_margin[6];
  /** offset between mesh lines of the last dimension */
  int q_2_off;
  /** offset between mesh lines of the two last dimensions */
  int q_21_off;
} local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  /** dimension of sub meshs to send. */
  int s_dim[6][3];
  /** left down corners of sub meshs to send. */
  int s_ld[6][3];
  /** up right corners of sub meshs to send. */
  int s_ur[6][3];
  /** sizes for send buffers. */
  int s_size[6];
  /** dimensionof sub meshs to recv. */
  int r_dim[6][3];
  /** left down corners of sub meshs to recv. */
  int r_ld[6][3];
  /** up right corners of sub meshs to recv. */
  int r_ur[6][3];
  /** sizes for recv buffers. */
  int r_size[6];
  /** maximal size for send/recv buffers. */
  int max;
} send_mesh;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** P3M parameters. */
extern p3m_struct p3m;

//The following are defined in p3m.c :

#define CA_INCREMENT 32  

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/// print the p3m parameters to the interpreters result
int printChargeP3MToResult(Tcl_Interp *interp);

/// print the p3m parameters to the interpreters result
int printDipolarP3MToResult(Tcl_Interp *interp);

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm.
 */
void P3M_init();

/** Calculate the k-space contribution to the coulomb interaction
     forces for electro- and/or magnetostatic interactions, if
    compiled in. */ 
double P3M_calc_kspace_forces(int force_flag, int energy_flag);

/** Clean up P3M memory allocations. */
void P3M_exit();

/*@}*/


/**************** headers for each specific p3m algorithm ********************************************/

#ifdef ELECTROSTATICS
#include "p3m-charges.h"  /* charge-charge interaction */
#endif

#ifdef MAGNETOSTATICS
#include "p3m-dipoles.h"  /* magnetic dipole-dipole interaction */
#endif 

#undef P3M_H_CURRENT

#endif /* of ifdef ELP3M */

#endif  /*of ifndef P3M_H */
