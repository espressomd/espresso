#ifndef P3M_H 
#define P3M_H 
/** \file p3m.h
    For more information about the p3m algorithm,
    see \ref p3m.c "p3m.c"
 */

/************************************************
 * data types
 ************************************************/

/** Structure to hold P3M parameters and some dependend variables. */
typedef struct {
  double bjerrum;   /** Bjerrum-length. */
  double alpha;     /** Ewald splitting parameter. */
  double r_cut;     /** Cutoff radius for real space electrostatics. */

  int    mesh[3];   /** number of mesh points per coordinate direction. */
  double mesh_off[3]; /** offset of the first mesh point (lower left 
			  corner) from the coordinate origin. */
  int    cao;       /** charge assignment order. */

  double epsilon;   /** epsilon of the "surrounding dielectric". */
  double prefactor; /** Coulomb prefactor, e.g. Bjerrum*Temp. */
  double r_cut2;    /** Cutoff radius squared. */
} p3m_struct;

/************************************************
 * exported variables
 ************************************************/

extern p3m_struct p3m;

/************************************************
 * exported functions
 ************************************************/

void   P3M_init();
void   P3M_perform();
void   P3M_exit();

/** Callback for setmd bjerrum. If the Bjerrum length is 0, the cutoff
    is also set to 0 to avoid unnecessary computation of the
    electrostatics. */
int bjerrum_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_alpha*/
int p3malpha_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_r_cut*/
int p3mrcut_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_mesh*/
int p3mmesh_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_cao*/
int p3mcao_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_epsilo*/
int p3mepsilon_callback(Tcl_Interp *interp, void *_data);
/** Callback for setmd p3m_mesh_offset*/
int p3mmeshoff_callback(Tcl_Interp *interp, void *_data);

#endif
