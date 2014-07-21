#ifndef VVOLUME_H
#define VVOLUME_H

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "integrate.hpp"

extern int setvo;
extern int vescnum;
extern double *CentVV;
extern double *VVol;
extern double VVolo[200];

#define ANGLE id1 == 1 && id2 == 0 && id3 == 3

//Calculate Centroid of Body/Vesicle
inline void GetCentroidV() {
    double temp[4*vescnum], upos[3];
    int img[3];
    int i, c, np, mid, midme;
    Cell *cell ;
    Particle *p ;
    
    for(i=0; i<4*vescnum; i++) {
	temp[i] = 0.0;
    }
    
    //All particles stored on local core
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {
			mid = p[i].p.mol_id;
	
	//if particle is part of a vesicle -> mol_id > 0
			if(mid > 0) {
	    
	    		midme = mid-1;
	    
	    		//Unfold Position
	    		upos[0] = p[i].r.p[0]; upos[1] = p[i].r.p[1]; upos[2] = p[i].r.p[2];
	    		img[0] = p[i].l.i[0]; img[1] = p[i].l.i[1]; img[2] = p[i].l.i[2];
	    
	    		unfold_position(upos,img);
	    
	    
	    //Add up all positions, a add 1 to number of particles
	    		temp[0+4*midme]+=upos[0];
	    		temp[1+4*midme]+=upos[1];
	    		temp[2+4*midme]+=upos[2];
	    		temp[3+4*midme]++;
	    		
			}
      }
    }
    
#ifdef VVOLUME_TRACE
    fprintf(stderr, "@node %d: t = %lf Before MPI_Allreduce Centroid\n", this_node, sim_time);
#endif    
    
    MPI_Allreduce(temp,CentVV,4*vescnum,MPI_DOUBLE,MPI_SUM,comm_cart);  
    
#ifdef VVOLUME_TRACE
    fprintf(stderr, "@node %d: t = %lf After MPI_Allreduce Centroid\n", this_node, sim_time);
#endif      
    
    //At this point I assume that there is indeed a body, otherwise we could divide by 0 if no particle is there
    for(i=0; i< vescnum; i++) {
    	CentVV[i*4]=CentVV[i*4]/CentVV[(i*4)+3];
    	CentVV[(i*4)+1]=CentVV[(i*4)+1]/CentVV[(i*4)+3];
    	CentVV[(i*4)+2]=CentVV[(i*4)+2]/CentVV[(i*4)+3];
    }
}

inline void GetVolumeV() {
    double temp[vescnum];
    double x1[3], x2[3], x3[3], rm[3], vprod[3];
    int i, c, np,j, type, n_partners, type_num, midme;
    double Vj;
    Cell *cell ;
    Particle *p ;
    Particle *p2, *p3=NULL, *p4=NULL;
    Bonded_ia_parameters *iaparams;
    char *errtxt;
    //int id1, id2, id3;
    
    for(i=0;i<vescnum;i++) {
    	temp[i]=0.0;
    }
    
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part;
      np = cell->n;

      for (i=0;i<np;i++) {
	//if particle is part of a vesicle -> mol_id > 0
	if(p[i].p.mol_id>0) {
	    j = 0;
	    midme = p[i].p.mol_id-1;
	    //Get Centroid
	    rm[0] = CentVV[midme*4];
	    rm[1] = CentVV[(midme*4)+1];
	    rm[2] = CentVV[(midme*4)+2];
	    
	    
	    
	    //Look at all bonds stored at that particle
	    while(j<p[i].bl.n) {
	      type_num = p[i].bl.e[j++];
	      iaparams = &bonded_ia_params[type_num];
	      type = iaparams->type;
	      n_partners = iaparams->num;
	      
	      //Partners needed anyway, since we need to go through bl.e
	      p2 = local_particles[p[i].bl.e[j++]];
	      if (!p2) {
	      errtxt = runtime_error(128 + 2*ES_INTEGER_SPACE);
	      ERROR_SPRINTF(errtxt,"{078 bond broken between particles %d and %d (particles not stored on the same node)} ",
	      p[i].p.identity, p[i].bl.e[i-1]);
	      return;
	      }

	      /* fetch particle 3 if needed */
	      if (n_partners >= 2) {
		p3 = local_particles[p[i].bl.e[j++]];
		if (!p3) {
		errtxt = runtime_error(128 + 3*ES_INTEGER_SPACE);
		ERROR_SPRINTF(errtxt,"{079 bond broken between particles %d, %d and %d (particles not stored on the same node)} ",
		p[i].p.identity, p[i].bl.e[j-2], p[i].bl.e[j-1]);
		return;
		}
	      }
	      
	      
	      /* fetch particle 4 if needed */
	      if (n_partners >= 3) {
	      p4 = local_particles[p[i].bl.e[j++]];
	      if (!p4) {
		errtxt = runtime_error(128 + 4*ES_INTEGER_SPACE);
		ERROR_SPRINTF(errtxt,"{080 bond broken between particles %d, %d, %d and %d (particles not stored on the same node)} ",
		p[i].p.identity, p[i].bl.e[j-3], p[i].bl.e[j-2], p[i].bl.e[j-1]);
		return;
		}
	      }
	      
	      //Get Volume Fragment
	      if(type == TRIEL_IA) {
		
		get_mi_vector(x1,p[i].r.p,rm);
		get_mi_vector(x2,p2->r.p,rm);
		get_mi_vector(x3,p3->r.p,rm);
		
		vector_product(x3,x2,vprod);
		Vj = -1.0*scalar(vprod,x1)/6.0;  
		temp[midme]+=Vj;
	      }
	      
	      
	    } //end loop over all bonds at partilce
	}//end if mol_id>0
      }//end for particles in cell
    } //end for all cells
     
#ifdef VVOLUME_TRACE
    fprintf(stderr, "@node %d: t = %lf Before MPI_Allreduce Volume\n", this_node, sim_time);
#endif    
     
     MPI_Allreduce(temp, VVol, vescnum, MPI_DOUBLE, MPI_SUM, comm_cart);

#ifdef VVOLUME_TRACE
    fprintf(stderr, "@node %d: t = %lf After MPI_Allreduce Volumne\n", this_node, sim_time);
#endif        

}

inline void RescaleVesicle() {
    int i, c, np, j;
    Cell *cell ;
    Particle *p ;
    double a;
    int midme;
    int img[3];
    double upos[3];
    double uposn[3], dp[3];
    
    double skin2 = SQR(0.5 * skin);
    
    
    //All particles stored on local core
    for (c=0;c<local_cells.n;c++) {
      cell = local_cells.cell[c] ;
      p = cell->part ;
      np = cell->n ;

      for (i=0;i<np;i++) {
	
	//if particle is part of a vesicle -> mol_id > 0
	if(p[i].p.mol_id>0) {
	    midme = p[i].p.mol_id-1;
	    a = pow(VVolo[midme]/VVol[midme], 1.0/3.0);
	    
	    //printf("@node %d: pid = %d, VVolo[0] = %lf, VVol[0] = %lf  a = %lf\n", this_node, p[i].p.identity, VVolo[0], VVol[0],a );
	    
	    //Unfold pos, since centroid in 'real' coordinates
	    upos[0] = p[i].r.p[0]; upos[1] = p[i].r.p[1]; upos[2] = p[i].r.p[2];
	    img[0] = p[i].l.i[0]; img[1] = p[i].l.i[1]; img[2] = p[i].l.i[2];
	    
	    unfold_position(upos,img);
	    
	    
	    //scale the position so that V = Vo again
	    for(j=0;j<3;j++) {
	     uposn[j] = a * upos[j] + (1-a) * CentVV[(midme*4)+j];
	     dp[j] = uposn[j]-upos[j];
	     p[i].r.p[j]+=dp[j];
	    }
	    
	    if(distance2(p[i].r.p,p[i].l.p_old) > skin2 ) resort_particles = 1;
	    
	 } //end if p of body
      } //end for all parts
    } // end for all cells
}


void SetCentVV();
void SetVVol();

inline void SetVo() {
	int i;
	for(i=0;i<vescnum;i++) {
		VVolo[i] = VVol[i];
	}
}

#endif
