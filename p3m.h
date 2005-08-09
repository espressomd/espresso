// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef P3M_H 
#define P3M_H
/** \file p3m.h   P3M algorithm for long range coulomb interaction.
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
 *  </ul>
 *
 *  For more information about the p3m algorithm,
 *  see \ref p3m.c "p3m.c"
 */

#include "utils.h"
#include "interaction_data.h"

#ifdef ELECTROSTATICS

/** This value for p3m.epsilon indicates metallic boundary conditions. */
#define P3M_EPSILON_METALLIC 0.0

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
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
} p3m_struct;

/** \name Exported Variables */
/************************************************************/
/*@{*/

/** P3M parameters. */
extern p3m_struct p3m;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/// print the p3m parameters to the interpreters result
int printP3MToResult(Tcl_Interp *interp);

/// parse the optimization parameters of p3m and the tuner
int inter_parse_p3m_opt_params(Tcl_Interp * interp, int argc, char ** argv);

/// parse the basic p3m parameters
int inter_parse_p3m(Tcl_Interp * interp, int argc, char ** argv);
  
/// sanity checks
int P3M_sanity_checks();

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm.
 */
void P3M_init();

/** Calculate number of charged particles, the sum of the squared
    charges and the squared sum of the charges. */
void P3M_count_charged_particles();

/** Updates \ref p3m_struct::alpha and \ref p3m_struct::r_cut if \ref box_l changed. */
void P3M_scaleby_box_l();

/** Tune P3M parameters to desired accuracy.

    Usage:
    \verbatim inter coulomb <bjerrum> p3m tune accuracy <value> [r_cut <value> mesh <value> cao <value>] \endverbatim

    The parameters are tuned to obtain the desired accuracy in best
    time, by running mpi_integrate(0) for several parameter sets.

    The function utilizes the analytic expression of the error estimate 
    for the P3M method in the book of Hockney and Eastwood (Eqn. 8.23) in 
    order to obtain the rms error in the force for a system of N randomly 
    distributed particles in a cubic box.
    For the real space error the estimate of Kolafa/Perram is used. 

    Parameter range if not given explicit values: For \ref p3m_struct::r_cut_iL
    the function uses the values (\ref min_local_box_l -\ref #skin) /
    (n * \ref box_l), n being an integer (this implies the assumption that \ref
    p3m_struct::r_cut_iL is the largest cutoff in the system!). For \ref
    p3m_struct::mesh the function uses the two values which matches best the
    equation: number of mesh point = number of charged particles. For
    \ref p3m_struct::cao the function considers all possible values.

    For each setting \ref p3m_struct::alpha_L is calculated assuming that the
    error contributions of real and reciprocal space should be equal.

    After checking if the total error fulfils the accuracy goal the
    time needed for one force calculation (including verlet list
    update) is measured via \ref mpi_integrate(0).

    The function returns a log of the performed tuning.

    The function is based on routines of the program HE_Q.c written by M. Deserno.
 */
int P3M_tune_parameters(Tcl_Interp *interp);

/** a probably faster adaptive tuning method. Uses the same error estimates and parameters as
    \ref P3M_adaptive_tune_parameters, but a different strategy for finding the optimum. The algorithm
    basically determines the mesh, cao and then the real space cutoff, in this nested order.

    For each mesh, the cao optimal for the mesh tested previously is used as an initial guess,
    and the algorithm tries whether increasing or decreasing it leads to a better solution. This
    is efficient, since the optimal cao only changes little with the meshes in general.

    The real space cutoff for a given mesh and cao is determined via a bisection on the error estimate,
    which determines where the error estimate equals the required accuracy. Therefore the smallest 
    possible, i.e. fastest real space cutoff is determined.

    Both the search over mesh and cao stop to search in a specific direction once the computation time is
    significantly higher than the currently known optimum.

    Compared to \ref P3M_tune_parameters, this function will test more parameters sets for efficiency, but
    the error estimate is calculated less often. In general this should be faster and give better results.
 */
int P3M_adaptive_tune_parameters(Tcl_Interp *interp);

/** Calculate the k-space contribution to the coulomb interaction
    forces. */ 
double P3M_calc_kspace_forces(int force_flag, int energy_flag);

/** Calculate real space contribution of coulomb pair forces.
    If NPT is compiled in, it returns the energy, which is needed for NPT. */
MDINLINE double add_p3m_coulomb_pair_force(Particle *p1, Particle *p2,
					   double *d,double dist2,double dist,double force[3])
{
  int j;
  double fac1,fac2, adist, erfc_part_ri;

  if(dist < p3m.r_cut) {
  if(dist > 0){		//Vincent
    adist = p3m.alpha * dist;
    erfc_part_ri = AS_erfc_part(adist) / dist;
    fac1 = coulomb.prefactor * p1->p.q * p2->p.q  * exp(-adist*adist);
    fac2 = fac1 * (erfc_part_ri + 2.0*p3m.alpha*wupii) / dist2;
    for(j=0;j<3;j++)
      force[j] += fac2 * d[j];
    ESR_TRACE(fprintf(stderr,"%d: RSE: Pair (%d-%d) dist=%.3f: force (%.3e,%.3e,%.3e)\n",this_node,
		      p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: ESR  f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: ESR  f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

#ifdef NPT
    return fac1 * erfc_part_ri;
#endif
  }}
  return 0.0;
}

/** Calculate real space contribution of coulomb pair energy. */
MDINLINE double p3m_coulomb_pair_energy(Particle *p1, Particle *p2,
				     double *d,double dist2,double dist)
{
  double adist, erfc_part_ri;

  if(dist < p3m.r_cut) {
    adist = p3m.alpha * dist;
    erfc_part_ri = AS_erfc_part(adist) / dist;
    return coulomb.prefactor*p1->p.q*p2->p.q *erfc_part_ri*exp(-adist*adist);
  }
  return 0.0;
}

#ifdef DIPOLES
/** Calculate real space contribution of dipolar pair energy. */
MDINLINE double p3m_dipol_pair_energy(Particle *p1, Particle *p2,
				     double *d,double dist2,double dist)
{
  double adist,adist2;
  double mu1_dot_mu2;
  double rij[3];
  double mu1_dot_rij;
  double mu2_dot_rij;
  int	i;
  double B,C;
  double erfc_ar;
  double dist3 = dist*dist2;
  


  if(dist < p3m.r_cut) {

    mu1_dot_mu2 = p1->r.dip[0]*p2->r.dip[0]
        	 +p1->r.dip[1]*p2->r.dip[1]
		 +p1->r.dip[2]*p2->r.dip[2];

    //Relative vector position:
    for (i=0;i<3;i++)
       rij[i] = p2->r.p[i] - p1->r.p[i];

    mu1_dot_rij =  p1->r.dip[0]*rij[0]
        	 + p1->r.dip[1]*rij[1]
		 + p1->r.dip[2]*rij[2];

    mu2_dot_rij =  p2->r.dip[0]*rij[0]
        	 + p2->r.dip[1]*rij[1]
		 + p2->r.dip[2]*rij[2];

    adist = p3m.alpha * dist;
    adist2 = SQR(adist);
    erfc_ar = AS_erfc_part(adist)*exp(-adist2);
    B = erfc_ar/dist3 + exp(-adist2)*2*p3m.alpha/(dist2*wupi);
    C = (3*erfc_ar + (2*p3m.alpha*dist*wupii)*(3+2*SQR(p3m.alpha*dist))*exp(-adist2))
        /(dist3*dist2);
    return coulomb.prefactor*(p1->p.dipm*p2->p.dipm)*
              (mu1_dot_mu2*B - mu1_dot_rij*mu2_dot_rij*C);
  }
  return 0.0;
}
#endif

/** Clean up P3M memory allocations. */
void   P3M_exit();

/*@}*/
#endif

#endif
