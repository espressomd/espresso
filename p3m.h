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

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
  double bjerrum;    /** Bjerrum-length (>0). */
  double alpha;      /** Ewald splitting parameter (0<alpha<1). */
  double r_cut;      /** Cutoff radius for real space electrostatics (>0). */

  int    mesh[3];    /** number of mesh points per coordinate direction (>0). */
  double mesh_off[3];/** offset of the first mesh point (lower left 
			 corner) from the coordinate origin ([0,1[). */
  int    cao;        /** charge assignment order ([0,7]). */
  int    inter;      /** number of interpolation points for charge assignment function */

  double epsilon;    /** epsilon of the "surrounding dielectric". */
  double prefactor;  /** Coulomb prefactor, e.g. Bjerrum*Temp. */
  double r_cut2;     /** Cutoff radius squared. */
  double cao_cut[3]; /** Cutoff for charge assignment. */
  double a[3];       /** mesh constant. */
  double ai[3];      /** inverse mesh constant. */
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

/** Initialize all structures, parameters and arrays needed for the 
 *  P3M algorithm.
 */
void   P3M_init();
/** Calculate the k-space contribution to the coulomb interaction forces. */ 
void   P3M_calc_kspace_forces();
/** Clean up P3M memory allocations. */
void   P3M_exit();

/** Callback for setmd bjerrum (bjerrum >= 0). 
    Bjerrum acts as switch to turn off all long range computations for bjerrum=0. */
int bjerrum_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_alpha (0.0 <= alpha <=1.0). */
int p3malpha_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_r_cut (r_cut >=0). */
int p3mrcut_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_mesh (all components > 0).
    So far only cubic meshs are supported. */
int p3mmesh_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_cao (1 <= cao <= 7, inter > 0). */
int p3mcao_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_epsilon. */
int p3mepsilon_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_mesh_offset (all components between 0.0 and 1.0). */
int p3mmeshoff_callback(Tcl_Interp *interp, void *_data);

/*@}*/

#endif
