/** \file p3m.c
 *
 *  Calculation of long range part of the coulomb interaction.
 *
 *  We use here a P3M (Particle-Particle Particle-Mesh) method based
 *  on the Ewald summation. Details of the used method can be found in
 *  Hockney/Eastwood and Deserno/Holm. The file p3m contains only the
 *  Particle-Mesh part.
 *
 *  The routines used in p3m.c are based on the work of M. Deserno and C. Holm.
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "debug.h"
#include "grid.h"
#include "integrate.h"
#include "particle_data.h"
#include "utils.h"
#include "communication.h"
#include "fft.h"
#include "p3m.h"

/************************************************
 * MACROS
 ************************************************/

/** number of Brillouin zones in the aliasing sums. */
#define BRILLOUIN 1       
/** increment size of charge assignment fields. */
#define CA_INCREMENT 10       

/************************************************
 * data types
 ************************************************/

/** Structure for local mesh parameters. */
typedef struct {
  /* local mesh characterization. */
  int dim[3];       /** dimension (size) of local mesh. */
  int size;         /** number of local mesh points. */
  int ld_ind[3];    /** index of lower left corner of the 
			local mesh in the global mesh. */
  double ld_pos[3]; /** position of the first local mesh point. */
  int inner[3];     /** dimension of mesh inside node domain. */
  int in_ld[3];     /** inner left down grid point */
  int in_ur[3];     /** inner up right grid point + (1,1,1)*/
  int margin[6];    /** number of margin mesh points. */
  int r_margin[6];  /** number of margin mesh points from neighbour nodes */
} local_mesh;

/** Structure for send/recv meshs. */
typedef struct {
  int s_dim[6][3];   /** dimension of sub meshs to send. */
  int s_ld[6][3];    /** left down corners of sub meshs to send. */
  int s_ur[6][3];    /** up right corners of sub meshs to send. */
  int s_size[6];     /** sizes for send buffers. */
  int r_dim[6][3];   /** dimensionof sub meshs to recv. */
  int r_ld[6][3];    /** left down corners of sub meshs to recv. */
  int r_ur[6][3];    /** up right corners of sub meshs to recv. */
  int r_size[6];     /** sizes for recv buffers. */
  int max;           /** maximal size for send/recv buffers. */
} send_mesh;

/************************************************
 * variables
 ************************************************/

/** p3m parameters. */
p3m_struct p3m;
/** local mesh. */
static local_mesh lm;
/** send/recv mesh sizes */
static send_mesh  sm;

/** size of linear array for local CA/FFT mesh . */
int    ca_mesh_size;
/** real space mesh (local) for CA/FFT.*/
double *rs_mesh;
/** k space mesh (local) for k space calculation and FFT.*/
double *ks_mesh;

/** Field to store grid points to send */
double *send_grid; 
/** Field to store grid points to recv */
double *recv_grid;

/** interpolation of the charge assignment function. */
double *int_caf[7];
/** position shift for calc. of first assignment mesh point. */
double pos_shift;

/** help variable for calculation of aliasing sums */
double *meshift;
/** Spatial differential operator in k-space. We use an i*k differentiation. */
double *d_op;
/** Optimal influence function (k-space) */
double *g;

/** number of charged particles on the node. */
int ca_num;
/** Charge fractions for mesh assignment. */
double *ca_frac;
/** first mesh point for charge assignment. */
int *ca_fmp;

/************************************************
 * privat functions
 ************************************************/

/** Calculates properties of the local FFT mesh for the charge assignment process. */
void calc_local_ca_mesh();
/** print local mesh content. */
void p3m_print_local_mesh(local_mesh l);
/** Calculates the properties of the send/recv sub-meshes of the local FFT mesh. 
 *  In order to calculate the recv sub-meshes there is a communication of 
 *  the margins between neighbouring nodes. */ 
void calc_send_mesh();
/** print send mesh content. */
void p3m_print_send_mesh(send_mesh sm);

void gather_fft_grid();
void spread_force_grid();
/** realloc charge assignment fields. */
void realloc_ca_fields(int newsize);
/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  ouput array with dimension dim[3] at start position start[3].  
*/
void add_block(double *in, double *out, int start[3], int size[3], int dim[3]);

/** Interpolates the P-th order charge assignment function from
 * Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
 * are also tabulated in Deserno/Holm. */
void interpolate_charge_assignment_function();

/** shifts the mesh points by mesh/2 */
void calc_meshift();

/** Calculates the Fourier transformed differential operator.  
 *  Remark: This is done on the level of n-vectors and not k-vectors,
 *           i.e. the prefactor i*2*PI/L is missing! */
void calc_differential_operator();
/** Calculates the optimal influence function of Hockney and Eastwood. 
 *
 *  Each node calculates only the values for its domain in k-space
 *  (see fft_plan[2].mesh and fft_plan[2].start).

 *  See also: Hockney/Eastwood 8-22 (p275). Note the somewhat
 *  different convention for the prefactors, which is described in
 *  Deserno/Holm. */
void calc_influence_function();

/** Calculates the aliasing sums for the optimal influence function.
 *
 * Calculates the aliasing sums in the nominator and denominator of
 * the expression for the optimal influence function (see
 * Hockney/Eastwood: 8-22, p. 275).  
 *
 * @param  n           n-vector for which the aliasing sum is to be performed.
 * @param  nominator   aliasing sums in the nominator.
 * @return denominator aliasing sum in the denominator
 */
MDINLINE double perform_aliasing_sums(int n[3], double nominator[3]);

/************************************************
 * public functions
 ************************************************/

void   P3M_init()
{
  int i,n;
  double temperature=1.0;

  if(p3m.bjerrum == 0.0) {       
    p3m.r_cut  = 0.0;
    p3m.r_cut2 = 0.0;
    if(this_node==0) 
      P3M_TRACE(fprintf(stderr,"0: P3M_init: Bjerrum length is zero.\n");
		fprintf(stderr,"   Electrostatics switched off!\n"));
  }
  else {  
    P3M_TRACE(fprintf(stderr,"%d: P3M_init: \n",this_node));
    /* parameter checks */
    if( (p3m.mesh[0] != p3m.mesh[1]) || (p3m.mesh[1] != p3m.mesh[2]) ) {
      if(this_node==0) {
	fprintf(stderr,"0: P3M_init: WARNING: Non cubic mesh found!\n");
	fprintf(stderr,"   use p3m.mesh[0]=%d as cubic mesh size!\n)",p3m.mesh[0]);
      }
      p3m.mesh[1] = p3m.mesh[0];
      p3m.mesh[2] = p3m.mesh[0];
    }
    if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
      if(this_node==0) {
	fprintf(stderr,"0: P3M_init: SERIOUS WARNING:\n"); 
	fprintf(stderr,"   No long range interactions for non cubic box.\n"); 
	fprintf(stderr,"   Switch off long range interactions! \n");
      }
      p3m.bjerrum =  0.0;
      p3m.r_cut  = 0.0;
      p3m.r_cut2 = 0.0;
      return ;
    }

    p3m.prefactor = p3m.bjerrum*temperature; 
    p3m.r_cut2    = p3m.r_cut*p3m.r_cut;
    for(i=0;i<3;i++) {
      p3m.ai[i]      = (double)p3m.mesh[i]/box_l[i]; 
      p3m.a[i]       = 1.0/p3m.ai[i];
      p3m.cao_cut[i] = p3m.a[i]*p3m.cao/2.0;
    }
    ca_frac = malloc(p3m.cao*p3m.cao*p3m.cao*CA_INCREMENT*sizeof(double));
    ca_fmp  = malloc(3*CA_INCREMENT*sizeof(int));

    calc_local_ca_mesh();
    P3M_TRACE(p3m_print_local_mesh(lm));
    calc_send_mesh();
    // DEBUG
    for(n=0;n<n_nodes;n++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if(n==this_node) P3M_TRACE(p3m_print_send_mesh(sm));
    }
    send_grid = malloc(sizeof(double)*sm.max);
    recv_grid = malloc(sizeof(double)*sm.max);

    interpolate_charge_assignment_function();
    /* position offset for calc. of fisrt meshpoint */
    pos_shift = (double)((p3m.cao-1)/2) - (p3m.cao%2)/2.0;
    P3M_TRACE(fprintf(stderr,"%d: pos_shift = %f\n",this_node,pos_shift)); 
 
    /* FFT */
    ca_mesh_size = fft_init(rs_mesh,lm.dim,lm.margin);
    rs_mesh = malloc(ca_mesh_size*sizeof(double));
    ks_mesh = malloc(ca_mesh_size*sizeof(double));
 
    /* k-space part: */
    calc_differential_operator();
    calc_influence_function();

    MPI_Barrier(MPI_COMM_WORLD);   
    P3M_TRACE(fprintf(stderr,"%d: p3m initialized\n",this_node));
  }
}

void   P3M_perform()
{
  int i,n,d,i0,i1,i2,ind,j[3];
  double q;
  /* position of a particle in local mesh units */
  double pos[3];
  /* index of first assignment mesh point; argument for interpolated caf */
  int first[3],arg[3];
  /* index, index jumps for rs_mesh array */
  int q_ind, q_m_off, q_s_off;
  /* (2*p3m.inter) + 1 */
  int inter2;
  /* tmp variables */
  double tmp0,tmp1;
  /* Energy */
  double energy=0.0;
  /* charged particle counter, charge fraction counter */
  int cp_cnt=0, cf_cnt=0;
  /* Prefactor for force */
  double force_prefac;

  MPI_Barrier(MPI_COMM_WORLD);   
  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: \n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);   

  /* prepare local FFT mesh */
  for(n=0; n<lm.size; n++) rs_mesh[n] = 0.0;
  q_m_off = (lm.dim[2] - p3m.cao);
  q_s_off = lm.dim[2] * (lm.dim[1] - p3m.cao);
  
  /* === charge assignment === */
  force_prefac = p3m.prefactor / (double)(p3m.mesh[0]*p3m.mesh[1]*p3m.mesh[2]);
  inter2 = (p3m.inter*2)+1;
  for(n=0; n<n_particles; n++) {
    if( (q=particles[n].q) != 0.0 ) {

      /* particle position in mesh coordinates */
      for(d=0;d<3;d++) {
	pos[d]   = (particles[n].p[d]-(lm.ld_pos[d]+pos_shift))*p3m.ai[d];
	first[d] = (int) pos[d];
	ca_fmp[(3*cp_cnt)+d] = first[d];
	arg[d]   = (int) ((pos[d]-first[d])*inter2);

	P3M_TRACE(if( pos[d]<0.0 ) 
		  fprintf(stderr,"%d: rs_mesh underflow! (P%d at %f)\n",
			  this_node,particles[n].identity,particles[n].p[d]));
	P3M_TRACE(if( (first[d]+p3m.cao-1) > lm.dim[d] )
		  fprintf(stderr,"%d: rs_mesh overflow!  (P%d at %f)\n",
			  this_node,particles[n].identity,particles[n].p[d]));
      }

      /* charge assignment */
      q_ind = first[2] + lm.dim[2]*(first[1] + (lm.dim[1]*first[0]));
      for(i0=0; i0<p3m.cao; i0++) {
	tmp0 = q * int_caf[i0][arg[0]];
	for(i1=0; i1<p3m.cao; i1++) {
	  tmp1 = tmp0 * int_caf[i1][arg[1]];
	  for(i2=0; i2<p3m.cao; i2++) {
	    ca_frac[cf_cnt] = tmp1 * int_caf[i2][arg[2]];
	    rs_mesh[q_ind++] += ca_frac[cf_cnt];
	    cf_cnt++;
	  }
	  q_ind += q_m_off;
	}
	q_ind += q_s_off;
      }
      cp_cnt++;
      if( (cp_cnt+1)>ca_num ) realloc_ca_fields(cp_cnt+1);
    }
  }
  if( (cp_cnt+CA_INCREMENT)<ca_num ) realloc_ca_fields(cp_cnt+CA_INCREMENT);

  /* Gather information for FFT grid inside the nodes domain (inner local mesh) */
  gather_fft_grid();

  /* === Perform forward 3D FFT (Charge Assignment Mesh) === */
  fft_perform_forw(rs_mesh);

  /* === K Space Calculations === */
  MPI_Barrier(MPI_COMM_WORLD);   
  P3M_TRACE(fprintf(stderr,"%d: p3m_perform: k-Space\n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);   

  /* Energy */
  ind=0;
  for(i=0; i<fft_plan[2].new_size; i++) 
    energy += g[i] * (SQR(rs_mesh[ind++])+SQR(rs_mesh[ind++]));
  /* Force preparation */
  ind = 0;
  for(i=0; i<fft_plan[2].new_size; i++) {
    ks_mesh[ind] = g[i] * rs_mesh[ind]; ind++;
    ks_mesh[ind] = g[i] * rs_mesh[ind]; ind++;
  }

  /* === 3 Fold backward 3D FFT (Force Component Meshs) === */

  /* Force component loop */
  for(d=0;d<3;d++) {  
    /* srqt(-1)*k differentiation */
    ind=0;
    for(j[0]=0; j[0]<fft_plan[2].new_mesh[0]; j[0]++) {
      for(j[1]=0; j[1]<fft_plan[2].new_mesh[1]; j[1]++) {
	for(j[2]=0; j[2]<fft_plan[2].new_mesh[2]; j[2]++) {
	  /* i*k*(Re+i*Im) = - Im*k + i*Re*k     (i=sqrt(-1)) */ 
	  rs_mesh[ind] = -(ks_mesh[ind+1]*d_op[ j[d]+fft_plan[2].start[d] ]); ind++;
	  rs_mesh[ind] =   ks_mesh[ind-1]*d_op[ j[d]+fft_plan[2].start[d] ];  ind++;
	}
      }
    }
    /* Back FFT force componenet mesh */
    fft_perform_back(rs_mesh);
    /* redistribute force componenet mesh */
    spread_force_grid();
    /* Assign force component from mesh to particle */
    cp_cnt=0; cf_cnt=0;
    for(n=0; n<n_particles; n++) { 
      if( (q=particles[n].q) != 0.0 ) {
	q_ind = ca_fmp[(3*cp_cnt)+2] + lm.dim[2]*(ca_fmp[(3*cp_cnt)+1] + (lm.dim[1]*ca_fmp[(3*cp_cnt)+0]));
	for(i0=0; i0<p3m.cao; i0++) {
	  for(i1=0; i1<p3m.cao; i1++) {
	    for(i2=0; i2<p3m.cao; i2++) {
	      particles[n].f[d] -= force_prefac*ca_frac[cf_cnt]*rs_mesh[q_ind++]; 
	      cf_cnt++;
	    }
	    q_ind += q_m_off;
	  }
	  q_ind += q_s_off;
	}
	cp_cnt++;
      }
    }
  }


}

void   P3M_exit()
{
  int i;
  /* free memory */
  free(ca_frac);
  free(ca_fmp);
  free(send_grid);
  free(recv_grid);
  free(rs_mesh);
  free(ks_mesh); 
  free(rs_mesh); 
  for(i=0; i<p3m.cao; i++) free(int_caf[i]);
}

void calc_local_ca_mesh() {
  int i;
  int ind[3];

  /* inner left down grid point (global index) */
  for(i=0;i<3;i++) lm.in_ld[i] = (int)ceil(my_left[i]*p3m.ai[i]-p3m.mesh_off[i]);
  /* inner up right grid point (global index) */
  for(i=0;i<3;i++) lm.in_ur[i] = (int)floor(my_right[i]*p3m.ai[i]-p3m.mesh_off[i]);
  /* correct roundof errors at up right boundary */
  for(i=0;i<3;i++) if((my_right[i]*p3m.ai[i]-p3m.mesh_off[i])-lm.in_ur[i]==0.0) lm.in_ur[i]--;
  /* inner grid dimensions */
  for(i=0;i<3;i++) lm.inner[i] = lm.in_ur[i] - lm.in_ld[i] + 1;
  /* index of left down grid point in global mesh */
  for(i=0;i<3;i++) 
    lm.ld_ind[i]=(int)ceil((my_left[i]-p3m.cao_cut[i]-skin)*p3m.ai[i]-p3m.mesh_off[i]);
  /* spacial position of left down mesh point */
  for(i=0;i<3;i++) lm.ld_pos[i] = lm.ld_ind[i]*p3m.a[i] + p3m.mesh_off[i];
  /* left down margin */
  for(i=0;i<3;i++) lm.margin[i*2] = lm.in_ld[i]-lm.ld_ind[i];
  /* up right grid point */
  for(i=0;i<3;i++) ind[i]=(int)floor((my_right[i]+p3m.cao_cut[i]+skin)*p3m.ai[i]-p3m.mesh_off[i]);
  /* correct roundof errors at up right boundary */
  for(i=0;i<3;i++)
    if(((my_right[i]+p3m.cao_cut[i]+skin)*p3m.ai[i]-p3m.mesh_off[i])-ind[i]==0) ind[i]--;
  /* up right margin */
  for(i=0;i<3;i++) lm.margin[(i*2)+1] = ind[i] - lm.in_ur[i];

  /* grid dimension */
  lm.size=1; 
  for(i=0;i<3;i++) {lm.dim[i] = ind[i] - lm.ld_ind[i] + 1; lm.size*=lm.dim[i];}
  /* reduce inner grid indices from global to local */
  for(i=0;i<3;i++) lm.in_ld[i] = lm.margin[i*2];
  for(i=0;i<3;i++) lm.in_ur[i] = lm.margin[i*2]+lm.inner[i];
}

void calc_send_mesh()
{
  int i,j,evenodd;
  int done[3]={0,0,0};
  MPI_Status status;
  /* send grids */
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      /* left */
      sm.s_ld[i*2][j] = 0 + done[j]*lm.margin[j*2];
      if(j==i) sm.s_ur[i*2][j] = lm.margin[j*2]; 
      else     sm.s_ur[i*2][j] = lm.dim[j]-done[j]*lm.margin[(j*2)+1];
      /* right */
      if(j==i) sm.s_ld[(i*2)+1][j] = lm.in_ur[j];
      else     sm.s_ld[(i*2)+1][j] = 0 + done[j]*lm.margin[j*2];
      sm.s_ur[(i*2)+1][j] = lm.dim[j] - done[j]*lm.margin[(j*2)+1];
    }   
    done[i]=1;
  }
  sm.max=0;
  for(i=0;i<6;i++) {
    sm.s_size[i] = 1;
    for(j=0;j<3;j++) {
      sm.s_dim[i][j] = sm.s_ur[i][j]-sm.s_ld[i][j];
      sm.s_size[i] *= sm.s_dim[i][j];
    }
    if(sm.s_size[i]>sm.max) sm.max=sm.s_size[i];
  }
  /* communication */
  for(i=0;i<6;i++) {
    if(i%2==0) j = i+1;
    else       j = i-1;
    if(node_neighbors[i] != this_node) {
      /* two step communication: first all even positions than all odd */
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[i/2]+evenodd)%2==0)
	  MPI_Send(&(lm.margin[i]), 1, MPI_INT, 
		   node_neighbors[i],0,MPI_COMM_WORLD);
	else
	  MPI_Recv(&(lm.r_margin[j]), 1, MPI_INT,
		   node_neighbors[j],0,MPI_COMM_WORLD,&status);    
      }
    }
    else {
      lm.r_margin[j] = lm.margin[i];
    }
  }
  /* recv grids */
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) {
      if(j==i) {
	sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j] + lm.margin[2*j];
	sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j] + lm.r_margin[2*j];
	sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j] - lm.r_margin[(2*j)+1];
	sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j] - lm.margin[(2*j)+1];
      }
      else {
	sm.r_ld[ i*2   ][j] = sm.s_ld[ i*2   ][j];
	sm.r_ur[ i*2   ][j] = sm.s_ur[ i*2   ][j];
	sm.r_ld[(i*2)+1][j] = sm.s_ld[(i*2)+1][j];
	sm.r_ur[(i*2)+1][j] = sm.s_ur[(i*2)+1][j];
      }
    }
  for(i=0;i<6;i++) {
    sm.r_size[i] = 1;
    for(j=0;j<3;j++) {
      sm.r_dim[i][j] = sm.r_ur[i][j]-sm.r_ld[i][j];
      sm.r_size[i] *= sm.r_dim[i][j];
    }
    if(sm.r_size[i]>sm.max) sm.max=sm.r_size[i];
  }
}

void gather_fft_grid()
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;

  MPI_Barrier(MPI_COMM_WORLD);   
  P3M_TRACE(fprintf(stderr,"%d: gather_fft_grid:\n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);   

  /* direction loop */
  for(s_dir=0; s_dir<6; s_dir++) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(sm.s_size[s_dir]>0) 
      pack_block(rs_mesh, send_grid, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1);
    /* communication */
    if(node_neighbors[s_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[s_dir/2]+evenodd)%2==0) {
	  if(sm.s_size[s_dir]>0) 
	    MPI_Send(send_grid, sm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], 0, MPI_COMM_WORLD);
	}
	else {
	  if(sm.r_size[r_dir]>0) 
	    MPI_Recv(recv_grid, sm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], 0, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = recv_grid;
      recv_grid = send_grid;
      send_grid = tmp_ptr;
    }
    /* add recv block */
    if(sm.r_size[r_dir]>0) {
      add_block(recv_grid, rs_mesh, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim); 
    }
  }
}

void spread_force_grid()
{
  int s_dir,r_dir,evenodd;
  MPI_Status status;
  double *tmp_ptr;
  MPI_Barrier(MPI_COMM_WORLD);   
  P3M_TRACE(fprintf(stderr,"%d: spread_force_grid:\n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);   

  /* direction loop */
  for(s_dir=5; s_dir>=0; s_dir--) {
    if(s_dir%2==0) r_dir = s_dir+1;
    else           r_dir = s_dir-1;
    /* pack send block */ 
    if(sm.s_size[s_dir]>0) 
      pack_block(rs_mesh, send_grid, sm.r_ld[r_dir], sm.r_dim[r_dir], lm.dim, 1);
    /* communication */
    if(node_neighbors[r_dir] != this_node) {
      for(evenodd=0; evenodd<2;evenodd++) {
	if((node_pos[r_dir/2]+evenodd)%2==0) {
	  if(sm.r_size[r_dir]>0) 
	    MPI_Send(send_grid, sm.r_size[r_dir], MPI_DOUBLE, 
		     node_neighbors[r_dir], 0, MPI_COMM_WORLD);
   	}
	else {
	  if(sm.s_size[s_dir]>0) 
	    MPI_Recv(recv_grid, sm.s_size[s_dir], MPI_DOUBLE, 
		     node_neighbors[s_dir], 0, MPI_COMM_WORLD, &status); 	    
	}
      }
    }
    else {
      tmp_ptr = recv_grid;
      recv_grid = send_grid;
      send_grid = tmp_ptr;
    }
    /* un pack recv block */
    if(sm.s_size[s_dir]>0) {
      unpack_block(recv_grid, rs_mesh, sm.s_ld[s_dir], sm.s_dim[s_dir], lm.dim, 1); 
    }
  }
}


void interpolate_charge_assignment_function()
{
  /* REMARK: This function is taken unchanged from polysim9 by M. Deserno. */
  double dInterpol=(double)p3m.inter, x;
  long   i;

  P3M_TRACE(fprintf(stderr,"%d - interpolating (%d) the order-%d charge assignment function\n",
		    this_node,p3m.inter,p3m.cao));

  for(i=0;i<p3m.cao;i++) 
    int_caf[i] = malloc(sizeof(double)*(2*p3m.inter+1));

  switch (p3m.cao) {
  case 1 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 1.0;
    }
  } break;
  case 2 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 0.5-x;
      int_caf[1][i+p3m.inter] = 0.5+x;
    }
  } break;
  case 3 : { 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = 0.5*SQR(0.5 - x);
      int_caf[1][i+p3m.inter] = 0.75 - SQR(x);
      int_caf[2][i+p3m.inter] = 0.5*SQR(0.5 + x);
    }
  } break;
  case 4 :{ 
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
      int_caf[1][i+p3m.inter] = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
      int_caf[2][i+p3m.inter] = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
      int_caf[3][i+p3m.inter] = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    }
  } break;
  case 5 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
      int_caf[1][i+p3m.inter] = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
      int_caf[2][i+p3m.inter] = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
      int_caf[3][i+p3m.inter] = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
      int_caf[4][i+p3m.inter] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    }
  } break;
  case 6 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
      int_caf[1][i+p3m.inter] = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
      int_caf[2][i+p3m.inter] = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
      int_caf[3][i+p3m.inter] = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
      int_caf[4][i+p3m.inter] = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
      int_caf[5][i+p3m.inter] = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    }
  } break;
  case 7 : {
    for (i=-p3m.inter; i<=p3m.inter; i++) {
      x=i/(2.0*dInterpol);
      int_caf[0][i+p3m.inter] = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
      int_caf[1][i+p3m.inter] = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
      int_caf[2][i+p3m.inter] = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
      int_caf[3][i+p3m.inter] = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
      int_caf[4][i+p3m.inter] = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
      int_caf[5][i+p3m.inter] = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
      int_caf[6][i+p3m.inter] = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    }
  } break;
  default :{
    fprintf(stderr,"%d: Charge assignment order %d unknown.\n",this_node,p3m.cao);
  }
  }
}

void calc_meshift(void)
{
  int i;
  double dmesh;

  dmesh = (double)p3m.mesh[0];
  meshift = malloc(p3m.mesh[0]*sizeof(double));

  for (i=0; i<p3m.mesh[0]; i++) meshift[i] = i - dround(i/dmesh)*dmesh; 
}

void calc_differential_operator()
{
  int i;
  double dmesh;

  dmesh = (double)p3m.mesh[0];
  d_op = malloc(p3m.mesh[0]*sizeof(double));

  for (i=0; i<p3m.mesh[0]; i++) 
    d_op[i] = (double)i - dround((double)i/dmesh)*dmesh;

  d_op[p3m.mesh[0]/2] = 0;
}

void calc_influence_function()
{
  int i,n[3],ind;
  int end[3];
  int size=1;
  double fak1,fak2;
  double nominator[3]={0.0,0.0,0.0},denominator=0.0;

  calc_meshift();

  for(i=0;i<3;i++) {
    size *= fft_plan[2].new_mesh[i];
    end[i] = fft_plan[2].start[i] + fft_plan[2].new_mesh[i];
  }
  g = malloc(size*sizeof(double));

  fak1  = p3m.mesh[0]*p3m.mesh[0]*p3m.mesh[0]*2.0/(box_l[0]*box_l[0]);

  for(n[0]=fft_plan[2].start[0]; n[0]<end[0]; n[0]++) 
    for(n[1]=fft_plan[2].start[1]; n[1]<end[1]; n[1]++) 
      for(n[2]=fft_plan[2].start[2]; n[2]<end[2]; n[2]++) {
	ind = (n[2]-fft_plan[2].start[2]) + fft_plan[2].new_mesh[2] * ((n[1]-fft_plan[2].start[1]) + (fft_plan[2].new_mesh[1]*(n[0]-fft_plan[2].start[0])));
	if( (n[0]==0) && (n[1]==0) && (n[2]==0) )
	  g[ind] = 0.0;
	else if( (n[0]%(p3m.mesh[0]/2)==0) && 
		 (n[1]%(p3m.mesh[0]/2)==0) && 
		 (n[2]%(p3m.mesh[0]/2)==0) )
	  g[ind] = 0.0;
	else {
	  denominator = perform_aliasing_sums(n,nominator);
	  fak2 =  d_op[n[0]]*nominator[0] + d_op[n[1]]*nominator[1] + d_op[n[2]]*nominator[2];  
	  fak2 /= ( ( SQR(d_op[n[0]])+SQR(d_op[n[1]])+SQR(d_op[n[2]]) ) * SQR(denominator) );
	  g[ind] = fak1*fak2;
	}
      }
}

MDINLINE double perform_aliasing_sums(int n[3], double nominator[3])
{
  int i;
  double denominator=0.0;
  /* lots of temporary variables... */
  double sx,sy,sz,f1,f2,f3,mx,my,mz,nmx,nmy,nmz,nm2,expo;
  double limit = 30;

  for(i=0;i<3;i++) nominator[i]=0.0;
  f1 = 1.0/(double)p3m.mesh[0];
  f2 = SQR(PI/(p3m.alpha*box_l[0]));

  for(mx = -BRILLOUIN; mx <= BRILLOUIN; mx++) {
    nmx = meshift[n[0]] + p3m.mesh[0]*mx;
    sx  = pow(sinc(f1*nmx),2.0*p3m.cao);
    for(my = -BRILLOUIN; my <= BRILLOUIN; my++) {
      nmy = meshift[n[1]] + p3m.mesh[0]*my;
      sy  = sx*pow(sinc(f1*nmy),2.0*p3m.cao);
      for(mz = -BRILLOUIN; mz <= BRILLOUIN; mz++) {
	nmz = meshift[n[2]] + p3m.mesh[0]*mz;
	sz  = sy*pow(sinc(f1*nmz),2.0*p3m.cao);
	
	nm2          =  SQR(nmx)+SQR(nmy)+SQR(nmz);
	expo         =  f2*nm2;
	f3           =  (expo<limit) ? sz*exp(-expo)/nm2 : 0.0;

	nominator[0] += f3*nmx; 
	nominator[1] += f3*nmy; 
	nominator[2] += f3*nmz; 
	denominator  += sz;
      }
    }
  }
  return denominator;
}

void realloc_ca_fields(int newsize)
{
  int incr = 0;
  if( newsize > ca_num ) incr = (newsize - ca_num)/CA_INCREMENT +1;
  else if( newsize < ca_num ) incr = (newsize - ca_num)/CA_INCREMENT;
  if(incr != 0) {
    ca_num += incr;
    if(ca_num<CA_INCREMENT) ca_num = CA_INCREMENT;
    ca_frac = (double *)realloc(ca_frac, p3m.cao*p3m.cao*p3m.cao*ca_num*sizeof(double));
    ca_fmp  = (int *)realloc(ca_fmp, 3*ca_num*sizeof(int));
  }  
}

/** Add values of a 3d-grid input block (size[3]) to values of 3d-grid
 *  ouput array with dimension dim[3] at start position start[3].  
*/
void add_block(double *in, double *out, int start[3], int size[3], int dim[3])
{
  /* fast,mid and slow changing indices */
  int f,m,s;
  /* linear index of in grid, linear index of out grid */
  int li_in=0,li_out=0;
  /* offsets for indizes in output grid */
  int m_out_offset,s_out_offset;

  li_out = start[2] + ( dim[2]*( start[1] + (dim[1]*start[0]) ) );
  m_out_offset  = dim[2] - size[2];
  s_out_offset  = (dim[2] * (dim[1] - size[1]));

  for(s=0 ;s<size[0]; s++) {
    for(m=0; m<size[1]; m++) {
      for(f=0; f<size[2]; f++) {
	out[li_out++] += in[li_in++];
      }
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }
}

/************************************************
 * Callback functions 
 ************************************************/

int bjerrum_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "Bjerrum length must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.bjerrum = data;
  mpi_bcast_parameter(FIELD_BJERRUM);
  return (TCL_OK);
}

int p3malpha_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0.0 || data > 1.0) {
    Tcl_AppendResult(interp, "P3M alpha must be in interval [0.0,1.0]", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.alpha = data;
  mpi_bcast_parameter(FIELD_P3M_ALPHA);
  return (TCL_OK);
}

int p3mrcut_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  if (data < 0) {
    Tcl_AppendResult(interp, "P3M r_cut must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.r_cut = data;
  mpi_bcast_parameter(FIELD_P3M_RCUT);
  return (TCL_OK);
}

int p3mmesh_callback(Tcl_Interp *interp, void *_data)
{
  int i, *data = (int *)_data;
  if (data[0] < 0 || data[1] < 0 || data[2] < 0) {
    Tcl_AppendResult(interp, "mesh sizes must be positiv.", (char *) NULL);
    return (TCL_ERROR);
  }
  for(i=0;i<3;i++) p3m.mesh[i] = data[i];
  mpi_bcast_parameter(FIELD_P3M_MESH);
  return (TCL_OK);
}

int p3mcao_callback(Tcl_Interp *interp, void *_data)
{
  int *data = (int *)_data;
  if (data[0] < 1 || data[0] > 7 || data[1] < 0) {
    Tcl_AppendResult(interp, "P3M CAO must be in interval [1,7] and inter > 0.", (char *) NULL);
    return (TCL_ERROR);
  }
  p3m.cao = data[0];
  p3m.inter = data[1];
  mpi_bcast_parameter(FIELD_P3M_CAO);
  return (TCL_OK);
}

int p3mepsilon_callback(Tcl_Interp *interp, void *_data)
{
  double data = *(double *)_data;
  p3m.epsilon = data;
  mpi_bcast_parameter(FIELD_P3M_EPSILON);
  return (TCL_OK);
}

int p3mmeshoff_callback(Tcl_Interp *interp, void *_data)
{
  int i;
  double *data = (double *)_data;
  if (data[0] < 0.0 || data[0] >= 1.0 || data[1] < 0.0 || data[1] >= 1.0 || data[2] < 0.0 || data[2] >= 1.0 ) {
    Tcl_AppendResult(interp, "P3M mesh offsets must be in interval [0.0,1.0[", (char *) NULL);
    return (TCL_ERROR);
  }  for(i=0;i<3;i++) p3m.mesh_off[i] = data[i];
  mpi_bcast_parameter(FIELD_P3M_MESH_OFF);
  return (TCL_OK);
}

/************************************************
 * Debug functions printing p3m structures 
 ************************************************/

void p3m_print_local_mesh(local_mesh l) 
{
  fprintf(stderr,"%d: local_mesh: dim=(%d,%d,%d), size=%d\n",this_node,
	  l.dim[0],l.dim[1],l.dim[2],l.size);
  fprintf(stderr,"    ld_ind=(%d,%d,%d), ld_pos=(%f,%f,%f)\n",
	  l.ld_ind[0],l.ld_ind[1],l.ld_ind[2],
	  l.ld_pos[0],l.ld_pos[1],l.ld_pos[2]);
  fprintf(stderr,"    inner=(%d,%d,%d) [(%d,%d,%d)-(%d,%d,%d)]\n",
	  l.inner[0],l.inner[1],l.inner[2],
	  l.in_ld[0],l.in_ld[1],l.in_ld[2],
	  l.in_ur[0],l.in_ur[1],l.in_ur[2]);
  fprintf(stderr,"    margin = (%d,%d,%d,  %d,%d,%d)\n",
	  l.margin[0],l.margin[1],l.margin[2],l.margin[3],l.margin[4],l.margin[5]);
  fprintf(stderr,"    r_margin=(%d,%d,%d,  %d,%d,%d)\n",
	  l.r_margin[0],l.r_margin[1],l.r_margin[2],l.r_margin[3],l.r_margin[4],l.r_margin[5]);
}

void p3m_print_send_mesh(send_mesh sm) 
{
  int i;
  fprintf(stderr,"%d: send_mesh: max=%d\n",this_node,sm.max);
  for(i=0;i<6;i++) {
    fprintf(stderr,"  dir=%d: s_dim (%d,%d,%d)  s_ld (%d,%d,%d) s_ur (%d,%d,%d) s_size=%d\n",i,sm.s_dim[i][0],sm.s_dim[i][1],sm.s_dim[i][2],sm.s_ld[i][0],sm.s_ld[i][1],sm.s_ld[i][2],sm.s_ur[i][0],sm.s_ur[i][1],sm.s_ur[i][2],sm.s_size[i]);
    fprintf(stderr,"         r_dim (%d,%d,%d)  r_ld (%d,%d,%d) r_ur (%d,%d,%d) r_size=%d\n",sm.r_dim[i][0],sm.r_dim[i][1],sm.r_dim[i][2],sm.r_ld[i][0],sm.r_ld[i][1],sm.r_ld[i][2],sm.r_ur[i][0],sm.r_ur[i][1],sm.r_ur[i][2],sm.r_size[i]);
  }
}

void p3m_print_p3m_struct(p3m_struct ps) {
  fprintf(stderr,"%d: p3m_struct: \n",this_node);
  fprintf(stderr,"   bjerrum=%f, alpha=%f, r_cut=%f, r_cut2=%f\n",
	  ps.bjerrum,ps.alpha,ps.r_cut,ps.r_cut2);
  fprintf(stderr,"   mesh=(%d,%d,%d), mesh_off=(%.4f,%.4f,%.4f)\n",
	  ps.mesh[0],ps.mesh[1],ps.mesh[2],
	  ps.mesh_off[0],ps.mesh_off[1],ps.mesh_off[2]);
  fprintf(stderr,"   cao=%d, inter=%d, epsilon=%f, prefactor=%f\n",
	  ps.cao,ps.inter,ps.epsilon,ps.prefactor);
  fprintf(stderr,"   cao_cut=(%f,%f,%f)\n",
	  ps.cao_cut[0],ps.cao_cut[1],ps.cao_cut[2]);
  fprintf(stderr,"   a=(%f,%f,%f), ai=(%f,%f,%f)\n",
	  ps.a[0],ps.a[1],ps.a[2],ps.ai[0],ps.ai[1],ps.ai[2]);
}
