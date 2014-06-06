/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file fft.cpp
 *
 *  Routines, row decomposition, data structures and communication for the 3D-FFT. 
 *
*/


#include "fft-dipolar.hpp"

#ifdef DP3M

#include <fftw3.h>
/* our remapping of malloc interferes with fftw3's name mangling. */
void *fftw_malloc(size_t n);

#include <mpi.h>
#include "communication.hpp"
#include "grid.hpp"
#include "fft-common.hpp"

/************************************************
 * variables
 ************************************************/

fft_data_struct dfft;

/** communicate the grid data according to the given fft_forw_plan. 
 * \param plan communication plan (see \ref fft_forw_plan).
 * \param in   input mesh.
 * \param out  output mesh.
*/
void dfft_forw_grid_comm(fft_forw_plan plan, double *in, double *out);

/** communicate the grid data according to the given fft_forw_plan/fft_bakc_plan. 
 * \param plan_f communication plan (see \ref fft_forw_plan).
 * \param plan_b additional back plan (see \ref fft_back_plan).
 * \param in     input mesh.
 * \param out    output mesh.
*/
void dfft_back_grid_comm(fft_forw_plan plan_f, fft_back_plan plan_b, double *in, double *out);


void dfft_pre_init()
{
  fft_common_pre_init(&dfft);
}

int dfft_init(double **data, 
	      int *local_mesh_dim, int *local_mesh_margin, 
	      int* global_mesh_dim, double *global_mesh_off,
	      int *ks_pnum)
{
  int i,j;
  /* helpers */
  int mult[3];

  int n_grid[4][3]; /* The four node grids. */
  int my_pos[4][3]; /* The position of this_node in the node grids. */
  int *n_id[4];     /* linear node identity lists for the node grids. */
  int *n_pos[4];    /* positions of nodes in the node grids. */
  /* FFTW WISDOM stuff. */
  char wisdom_file_name[255];
  FILE *wisdom_file;
  int wisdom_status;

  FFT_TRACE(fprintf(stderr,"%d: dipolar dfft_init():\n",this_node));


  dfft.max_comm_size=0; dfft.max_mesh_size=0;
  for(i=0;i<4;i++) {
    n_id[i]  = (int *) malloc(1*n_nodes*sizeof(int));
    n_pos[i] = (int *) malloc(3*n_nodes*sizeof(int));
  }

  /* === node grids === */
  /* real space node grid (n_grid[0]) */
  for(i=0;i<3;i++) {
    n_grid[0][i] = node_grid[i];
    my_pos[0][i] = node_pos[i];
  }
  for(i=0;i<n_nodes;i++) {
    map_node_array(i,&(n_pos[0][3*i+0]));
    n_id[0][get_linear_index( n_pos[0][3*i+0],n_pos[0][3*i+1],n_pos[0][3*i+2], n_grid[0])] = i;
  }
    
  /* FFT node grids (n_grid[1 - 3]) */
  calc_2d_grid(n_nodes,n_grid[1]);
  /* resort n_grid[1] dimensions if necessary */
  dfft.plan[1].row_dir = map_3don2d_grid(n_grid[0], n_grid[1], mult);
  dfft.plan[0].n_permute = 0;
  for(i=1;i<4;i++) dfft.plan[i].n_permute = (dfft.plan[1].row_dir+i)%3;
  for(i=0;i<3;i++) {
    n_grid[2][i] = n_grid[1][(i+1)%3];
    n_grid[3][i] = n_grid[1][(i+2)%3];
  }
  dfft.plan[2].row_dir = (dfft.plan[1].row_dir-1)%3;
  dfft.plan[3].row_dir = (dfft.plan[1].row_dir-2)%3;



  /* === communication groups === */
  /* copy local mesh off real space charge assignment grid */
  for(i=0;i<3;i++) dfft.plan[0].new_mesh[i] = local_mesh_dim[i];
  for(i=1; i<4;i++) {
    dfft.plan[i].g_size=fft_find_comm_groups(n_grid[i-1], n_grid[i], n_id[i-1], n_id[i], 
					dfft.plan[i].group, n_pos[i], my_pos[i]);
    if(dfft.plan[i].g_size==-1) {
      /* try permutation */
      j = n_grid[i][(dfft.plan[i].row_dir+1)%3];
      n_grid[i][(dfft.plan[i].row_dir+1)%3] = n_grid[i][(dfft.plan[i].row_dir+2)%3];
      n_grid[i][(dfft.plan[i].row_dir+2)%3] = j;
      dfft.plan[i].g_size=fft_find_comm_groups(n_grid[i-1], n_grid[i], n_id[i-1], n_id[i], 
					  dfft.plan[i].group, n_pos[i], my_pos[i]);
      if(dfft.plan[i].g_size==-1) {
	fprintf(stderr,"%d: dipolar INTERNAL ERROR: fft_find_comm_groups error\n", this_node);
	errexit();
      }
    }

    dfft.plan[i].send_block = (int *)realloc(dfft.plan[i].send_block, 6*dfft.plan[i].g_size*sizeof(int));
    dfft.plan[i].send_size  = (int *)realloc(dfft.plan[i].send_size, 1*dfft.plan[i].g_size*sizeof(int));
    dfft.plan[i].recv_block = (int *)realloc(dfft.plan[i].recv_block, 6*dfft.plan[i].g_size*sizeof(int));
    dfft.plan[i].recv_size  = (int *)realloc(dfft.plan[i].recv_size, 1*dfft.plan[i].g_size*sizeof(int));

    dfft.plan[i].new_size = fft_calc_local_mesh(my_pos[i], n_grid[i], global_mesh_dim,
					   global_mesh_off, dfft.plan[i].new_mesh, 
					   dfft.plan[i].start);  
    permute_ifield(dfft.plan[i].new_mesh,3,-(dfft.plan[i].n_permute));
    permute_ifield(dfft.plan[i].start,3,-(dfft.plan[i].n_permute));
    dfft.plan[i].n_ffts = dfft.plan[i].new_mesh[0]*dfft.plan[i].new_mesh[1];

    /* === send/recv block specifications === */
    for(j=0; j<dfft.plan[i].g_size; j++) {
      int k, node;
      /* send block: this_node to comm-group-node i (identity: node) */
      node = dfft.plan[i].group[j];
      dfft.plan[i].send_size[j] 
	= fft_calc_send_block(my_pos[i-1], n_grid[i-1], &(n_pos[i][3*node]), n_grid[i],
			      global_mesh_dim, global_mesh_off, &(dfft.plan[i].send_block[6*j]));
      permute_ifield(&(dfft.plan[i].send_block[6*j]),3,-(dfft.plan[i-1].n_permute));
      permute_ifield(&(dfft.plan[i].send_block[6*j+3]),3,-(dfft.plan[i-1].n_permute));
      if(dfft.plan[i].send_size[j] > dfft.max_comm_size) 
	dfft.max_comm_size = dfft.plan[i].send_size[j];
      /* First plan send blocks have to be adjusted, since the CA grid
	 may have an additional margin outside the actual domain of the
	 node */
      if(i==1) {
	for(k=0;k<3;k++) 
	  dfft.plan[1].send_block[6*j+k  ] += local_mesh_margin[2*k];
      }
      /* recv block: this_node from comm-group-node i (identity: node) */
      dfft.plan[i].recv_size[j] 
	= fft_calc_send_block(my_pos[i], n_grid[i], &(n_pos[i-1][3*node]), n_grid[i-1],
			      global_mesh_dim, global_mesh_off,&(dfft.plan[i].recv_block[6*j]));
      permute_ifield(&(dfft.plan[i].recv_block[6*j]),3,-(dfft.plan[i].n_permute));
      permute_ifield(&(dfft.plan[i].recv_block[6*j+3]),3,-(dfft.plan[i].n_permute));
      if(dfft.plan[i].recv_size[j] > dfft.max_comm_size) 
	dfft.max_comm_size = dfft.plan[i].recv_size[j];
    }

    for(j=0;j<3;j++) dfft.plan[i].old_mesh[j] = dfft.plan[i-1].new_mesh[j];
    if(i==1) 
      dfft.plan[i].element = 1; 
    else {
      dfft.plan[i].element = 2;
      for(j=0; j<dfft.plan[i].g_size; j++) {
	dfft.plan[i].send_size[j] *= 2;
	dfft.plan[i].recv_size[j] *= 2;
      }
    }
    /* DEBUG */
    for(j=0;j<n_nodes;j++) {
      /* MPI_Barrier(comm_cart); */
      if(j==this_node) FFT_TRACE(fft_print_fft_plan(dfft.plan[i]));
    }
  }

  /* Factor 2 for complex fields */
  dfft.max_comm_size *= 2;
  dfft.max_mesh_size = (local_mesh_dim[0]*local_mesh_dim[1]*local_mesh_dim[2]);
  for(i=1;i<4;i++) 
    if(2*dfft.plan[i].new_size > dfft.max_mesh_size) dfft.max_mesh_size = 2*dfft.plan[i].new_size;

  FFT_TRACE(fprintf(stderr,"%d: dfft.max_comm_size = %d, dfft.max_mesh_size = %d\n",
		    this_node,dfft.max_comm_size,dfft.max_mesh_size));

  /* === pack function === */
  for(i=1;i<4;i++) {
    dfft.plan[i].pack_function = fft_pack_block_permute2; 
    FFT_TRACE(fprintf(stderr,"%d: forw plan[%d] permute 2 \n",this_node,i));
  }
  (*ks_pnum)=6;
  if(dfft.plan[1].row_dir==2) {
    dfft.plan[1].pack_function = fft_pack_block;
    FFT_TRACE(fprintf(stderr,"%d: forw plan[%d] permute 0 \n",this_node,1));
    (*ks_pnum)=4;
  }
  else if(dfft.plan[1].row_dir==1) {
    dfft.plan[1].pack_function = fft_pack_block_permute1;
    FFT_TRACE(fprintf(stderr,"%d: forw plan[%d] permute 1 \n",this_node,1));
    (*ks_pnum)=5;
  }
  
  /* Factor 2 for complex numbers */
  dfft.send_buf = (double *)realloc(dfft.send_buf, dfft.max_comm_size*sizeof(double));
  dfft.recv_buf = (double *)realloc(dfft.recv_buf, dfft.max_comm_size*sizeof(double));
  (*data)  = (double *)realloc((*data), dfft.max_mesh_size*sizeof(double));
  dfft.data_buf = (double *)realloc(dfft.data_buf, dfft.max_mesh_size*sizeof(double));
  if(!(*data) || !dfft.data_buf || !dfft.recv_buf || !dfft.send_buf) {
    fprintf(stderr,"%d: Could not allocate FFT data arays\n",this_node);
    errexit();
  }

  fftw_complex *c_data     = (fftw_complex *) (*data);

  /* === FFT Routines (Using FFTW / RFFTW package)=== */
  for(i=1;i<4;i++) {
    dfft.plan[i].dir = FFTW_FORWARD;   
    /* FFT plan creation. 
       Attention: destroys contents of c_data/data and c_data_buf/data_buf. */
    wisdom_status   = FFTW_FAILURE;
    sprintf(wisdom_file_name,"dfftw3_1d_wisdom_forw_n%d.file",
	    dfft.plan[i].new_mesh[2]);
    if( (wisdom_file=fopen(wisdom_file_name,"r"))!=NULL ) {
      wisdom_status = fftw_import_wisdom_from_file(wisdom_file);
      fclose(wisdom_file);
    }
    if(dfft.init_tag==1) fftw_destroy_plan(dfft.plan[i].our_fftw_plan);
//printf("dfft.plan[%d].n_ffts=%d\n",i,dfft.plan[i].n_ffts);
    dfft.plan[i].our_fftw_plan =
      fftw_plan_many_dft(1,&dfft.plan[i].new_mesh[2],dfft.plan[i].n_ffts,
                         c_data,NULL,1,dfft.plan[i].new_mesh[2],
                         c_data,NULL,1,dfft.plan[i].new_mesh[2],
                         dfft.plan[i].dir,FFTW_PATIENT);
    if( wisdom_status == FFTW_FAILURE && 
	(wisdom_file=fopen(wisdom_file_name,"w"))!=NULL ) {
      fftw_export_wisdom_to_file(wisdom_file);
      fclose(wisdom_file);
    }
    dfft.plan[i].fft_function = fftw_execute;        
  }

  /* === The BACK Direction === */
  /* this is needed because slightly different functions are used */
  for(i=1;i<4;i++) {
    dfft.back[i].dir = FFTW_BACKWARD;
    wisdom_status   = FFTW_FAILURE;
    sprintf(wisdom_file_name,"dfftw3_1d_wisdom_back_n%d.file",
	    dfft.plan[i].new_mesh[2]);
    if( (wisdom_file=fopen(wisdom_file_name,"r"))!=NULL ) {
      wisdom_status = fftw_import_wisdom_from_file(wisdom_file);
      fclose(wisdom_file);
    }    
    if(dfft.init_tag==1) fftw_destroy_plan(dfft.back[i].our_fftw_plan);
    dfft.back[i].our_fftw_plan =
      fftw_plan_many_dft(1,&dfft.plan[i].new_mesh[2],dfft.plan[i].n_ffts,
                         c_data,NULL,1,dfft.plan[i].new_mesh[2],
                         c_data,NULL,1,dfft.plan[i].new_mesh[2],
                         dfft.back[i].dir,FFTW_PATIENT);
    if( wisdom_status == FFTW_FAILURE && 
	(wisdom_file=fopen(wisdom_file_name,"w"))!=NULL ) {
      fftw_export_wisdom_to_file(wisdom_file);
      fclose(wisdom_file);
    }
    dfft.back[i].fft_function = fftw_execute;
    dfft.back[i].pack_function = fft_pack_block_permute1;
    FFT_TRACE(fprintf(stderr,"%d: back plan[%d] permute 1 \n",this_node,i));
  }
  if(dfft.plan[1].row_dir==2) {
    dfft.back[1].pack_function = fft_pack_block;
    FFT_TRACE(fprintf(stderr,"%d: back plan[%d] permute 0 \n",this_node,1));
  }
  else if(dfft.plan[1].row_dir==1) {
    dfft.back[1].pack_function = fft_pack_block_permute2;
    FFT_TRACE(fprintf(stderr,"%d: back plan[%d] permute 2 \n",this_node,1));
  }
  dfft.init_tag=1;
  /* free(data); */
  for(i=0;i<4;i++) { free(n_id[i]); free(n_pos[i]); }
  return dfft.max_mesh_size; 
}


void dfft_perform_forw(double *data)
{
  int i;
  /* int m,n,o; */
  /* ===== first direction  ===== */
  FFT_TRACE(fprintf(stderr,"%d: dipolar fft_perform_forw: dir 1:\n",this_node));

  fftw_complex *c_data     = (fftw_complex *) data;
  fftw_complex *c_data_buf = (fftw_complex *) dfft.data_buf;

  /* communication to current dir row format (in is data) */
  dfft_forw_grid_comm(dfft.plan[1], data, dfft.data_buf);


  /*
    fprintf(stderr,"%d: start grid \n",this_node);
    i=0;
    for(m=0;m<8;m++) {
    for(n=0;n<8;n++) {
    for(o=0;o<8;o++) {
    fprintf(stderr,"%.3f ",data_buf[i++]);
    }
    fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
    }
  */

  /* complexify the real data array (in is data_buf) */
  for(i=0;i<dfft.plan[1].new_size;i++) {
    data[2*i]     = dfft.data_buf[i];     /* real value */
    data[(2*i)+1] = 0;       /* complex value */
  }
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(dfft.plan[1].our_fftw_plan,c_data,c_data);
  /* ===== second direction ===== */
  FFT_TRACE(fprintf(stderr,"%d: dipolar fft_perform_forw: dir 2:\n",this_node));
  /* communication to current dir row format (in is data) */
  dfft_forw_grid_comm(dfft.plan[2], data, dfft.data_buf);
  /* perform FFT (in/out is data_buf)*/
  fftw_execute_dft(dfft.plan[2].our_fftw_plan,c_data_buf,c_data_buf);
  /* ===== third direction  ===== */
  FFT_TRACE(fprintf(stderr,"%d: dipolar fft_perform_forw: dir 3:\n",this_node));
  /* communication to current dir row format (in is data_buf) */
  dfft_forw_grid_comm(dfft.plan[3], dfft.data_buf, data);
  /* perform FFT (in/out is data)*/
  fftw_execute_dft(dfft.plan[3].our_fftw_plan,c_data,c_data);
  //fft_print_global_fft_mesh(dfft.plan[3],data,1,0);

  /* REMARK: Result has to be in data. */
}


void dfft_perform_back(double *data)
{
  int i;

  fftw_complex *c_data     = (fftw_complex *) data;
  fftw_complex *c_data_buf = (fftw_complex *) dfft.data_buf;
  
  /* ===== third direction  ===== */
  FFT_TRACE(fprintf(stderr,"%d: dipolar fft_perform_back: dir 3:\n",this_node));


  /* perform FFT (in is data) */
  fftw_execute_dft(dfft.back[3].our_fftw_plan,c_data,c_data);
  /* communicate (in is data)*/
  dfft_back_grid_comm(dfft.plan[3],dfft.back[3],data,dfft.data_buf);
 
  /* ===== second direction ===== */
  FFT_TRACE(fprintf(stderr,"%d: dipolar fft_perform_back: dir 2:\n",this_node));
  /* perform FFT (in is data_buf) */
  fftw_execute_dft(dfft.back[2].our_fftw_plan,c_data_buf,c_data_buf);
  /* communicate (in is data_buf) */
  dfft_back_grid_comm(dfft.plan[2],dfft.back[2],dfft.data_buf,data);

  /* ===== first direction  ===== */
  FFT_TRACE(fprintf(stderr,"%d: fft_perform_back: dir 1:\n",this_node));
  /* perform FFT (in is data) */
  fftw_execute_dft(dfft.back[1].our_fftw_plan,c_data,c_data);
  /* throw away the (hopefully) empty complex component (in is data)*/
  for(i=0;i<dfft.plan[1].new_size;i++) {
    dfft.data_buf[i] = data[2*i]; /* real value */
    //Vincent:
    if (data[2*i+1]>1e-5) {
      printf("dipoar fft - Complex value is not zero (i=%d,data=%g)!!!\n",i,data[2*i+1]);
      if (i>100) exit(-1);
    }
  }
  /* communicate (in is data_buf) */
  dfft_back_grid_comm(dfft.plan[1],dfft.back[1],dfft.data_buf,data);


  /* REMARK: Result has to be in data. */
}

void dfft_forw_grid_comm(fft_forw_plan plan, double *in, double *out)
{
  int i;
  MPI_Status status;
  double *tmp_ptr;

  for(i=0;i<plan.g_size;i++) {   
    plan.pack_function(in, dfft.send_buf, &(plan.send_block[6*i]), 
		       &(plan.send_block[6*i+3]), plan.old_mesh, plan.element);

    if(plan.group[i]<this_node) {       /* send first, receive second */
      MPI_Send(dfft.send_buf, plan.send_size[i], MPI_DOUBLE, 
	       plan.group[i], REQ_FFT_FORW, comm_cart);
      MPI_Recv(dfft.recv_buf, plan.recv_size[i], MPI_DOUBLE, 
	       plan.group[i], REQ_FFT_FORW, comm_cart, &status); 	
    }
    else if(plan.group[i]>this_node) {  /* receive first, send second */
      MPI_Recv(dfft.recv_buf, plan.recv_size[i], MPI_DOUBLE, 
	       plan.group[i], REQ_FFT_FORW, comm_cart, &status); 	
      MPI_Send(dfft.send_buf, plan.send_size[i], MPI_DOUBLE, 
	       plan.group[i], REQ_FFT_FORW, comm_cart);      
    }
    else {                              /* Self communication... */   
      tmp_ptr  = dfft.send_buf;
      dfft.send_buf = dfft.recv_buf;
      dfft.recv_buf = tmp_ptr;
    }
    fft_unpack_block(dfft.recv_buf, out, &(plan.recv_block[6*i]), 
		 &(plan.recv_block[6*i+3]), plan.new_mesh, plan.element);
  }
}

void dfft_back_grid_comm(fft_forw_plan plan_f,  fft_back_plan plan_b, double *in, double *out)
{
  int i;
  MPI_Status status;
  double *tmp_ptr;

  /* Back means: Use the send/recieve stuff from the forward plan but
     replace the recieve blocks by the send blocks and vice
     versa. Attention then also new_mesh and old_mesh are exchanged */

  for(i=0;i<plan_f.g_size;i++) {
    
    plan_b.pack_function(in, dfft.send_buf, &(plan_f.recv_block[6*i]), 
		       &(plan_f.recv_block[6*i+3]), plan_f.new_mesh, plan_f.element);

    if(plan_f.group[i]<this_node) {       /* send first, receive second */
      MPI_Send(dfft.send_buf, plan_f.recv_size[i], MPI_DOUBLE, 
	       plan_f.group[i], REQ_FFT_BACK, comm_cart);
      MPI_Recv(dfft.recv_buf, plan_f.send_size[i], MPI_DOUBLE, 
	       plan_f.group[i], REQ_FFT_BACK, comm_cart, &status); 	
    }
    else if(plan_f.group[i]>this_node) {  /* receive first, send second */
      MPI_Recv(dfft.recv_buf, plan_f.send_size[i], MPI_DOUBLE, 
	       plan_f.group[i], REQ_FFT_BACK, comm_cart, &status); 	
      MPI_Send(dfft.send_buf, plan_f.recv_size[i], MPI_DOUBLE, 
	       plan_f.group[i], REQ_FFT_BACK, comm_cart);      
    }
    else {                                /* Self communication... */   
      tmp_ptr  = dfft.send_buf;
      dfft.send_buf = dfft.recv_buf;
      dfft.recv_buf = tmp_ptr;
    }
    fft_unpack_block(dfft.recv_buf, out, &(plan_f.send_block[6*i]), 
		 &(plan_f.send_block[6*i+3]), plan_f.old_mesh, plan_f.element);
  }
}


#endif /* DP3M */
