#include "utils.h"

#ifdef COMFORCE
MDINLINE void calc_comforce()
{
  int t0,t1,k, j;
  IA_parameters *ia_params;
  double com0[3], com1[3], gyrtensor[9], diff[3];
  double vect0[3], vect1[3], eva[3], eve[3], fvect[3];
  Particle *p;
  int i, np, c;
  Cell *cell;
  
  for (t0=0; t0<n_particle_types-1; t0++) {
    for (t1=t0+1; t1<n_particle_types; t1++) {
      ia_params = get_ia_param(t0,t1);
      if(ia_params->COMFORCE_flag == 1) {
	      centermass(t0,com0);
	      centermass(t1,com1);
	      for (i = 0; i < 3; i++) {
		      diff[i]=com1[i]-com0[i];
	      }
        gyrationtensor(t0, gyrtensor);
        k=calc_eigenvalues_3x3(gyrtensor, eva);
        /* perpendicular force */
        if(ia_params->COMFORCE_dir == 1) {
	        k=calc_eigenvector_3x3(gyrtensor,eva[0],eve);
          /*By doing two vector products find radial axis along the target system */
	        vector_product(eve,diff,vect0);
	        vector_product(vect0,eve,vect1);
  
          /* normalize vect1, return is fvect */
	        unit_vector(vect1,fvect);
        } else {
        /* parallel force */
	        k=calc_eigenvector_3x3(gyrtensor,eva[0],fvect);
        }
        
        /* orient it along the com vector */
        if (scalar(fvect,diff) < 0.) {
	        for (i = 0; i < 3; i++) {
		        fvect[i] = -fvect[i];
	        }
        }
        
        /* Now apply the force */
        for (c = 0; c < local_cells.n; c++) {
          cell = local_cells.cell[c];
          p  = cell->part;
          np = cell->n;
          for(i = 0; i < np; i++) {
	 	        if(p[i].p.type==t0) {
      	      for(j = 0; j < 3; j++) {
				        p[i].f.f[j] -= ia_params->COMFORCE_fratio * ia_params->COMFORCE_force * fvect[j];
			        }
		        }
	 	        if(p[i].p.type==t1) {
      	      for(j = 0; j < 3; j++) {
				        p[i].f.f[j] +=  ia_params->COMFORCE_force * fvect[j];
			        }
		        }
          }
        }
        /*end of force application */
      }
    }
  }

}
#endif
