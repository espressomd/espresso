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
/** \file fft-common.cpp
 *
 *  Routines, row decomposition, data structures and communication for the 3D-FFT. 
 *
*/
#include "fft-common.hpp"

#if defined(P3M) || defined(DP3M)

#include <cstring>
#include <fftw3.h>
#include <mpi.h>

#include "utils.hpp"
#include "communication.hpp"


void fft_common_pre_init(fft_data_struct *fft)
{
  for(int i=0;i<4;i++) {
    fft->plan[i].group = (int*)malloc(1*n_nodes*sizeof(int));
    fft->plan[i].send_block = NULL;
    fft->plan[i].send_size  = NULL;
    fft->plan[i].recv_block = NULL;
    fft->plan[i].recv_size  = NULL;
  }

  fft->init_tag = 0;
  fft->max_comm_size = 0;
  fft->max_mesh_size = 0;
  fft->send_buf = NULL;
  fft->recv_buf = NULL;
  fft->data_buf = NULL;
}

void fft_pack_block(double *in, double *out, int start[3], int size[3], int dim[3], int element)
{
  /* mid and slow changing indices */
  int m,s;
  /* linear index of in grid, linear index of out grid */
  int li_in,li_out=0;
  /* copy size */
  int copy_size;
  /* offsets for indizes in input grid */
  int m_in_offset,s_in_offset;
  /* offsets for indizes in output grid */
  int m_out_offset;

  copy_size    = element * size[2] * sizeof(double);
  m_in_offset  = element * dim[2];
  s_in_offset  = element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = element * size[2];
  li_in        = element * (start[2]+dim[2]*(start[1]+dim[1]*start[0]));

  for(s=0 ;s<size[0]; s++) {
    for(m=0; m<size[1]; m++) {
      memcpy(&(out[li_out]), &(in[li_in]), copy_size);
      li_in  += m_in_offset;
      li_out += m_out_offset;
    }
    li_in += s_in_offset;
  }
}

void fft_pack_block_permute1(double *in, double *out, int start[3], int size[3], 
			 int dim[3], int element)
{
  /* slow,mid and fast changing indices for input  grid */
  int s,m,f,e;
  /* linear index of in grid, linear index of out grid */
  int li_in,li_out=0;
  /* offsets for indizes in input grid */
  int m_in_offset,s_in_offset;
  /* offset for mid changing indices of output grid */
  int m_out_offset;

  m_in_offset  =  element * (dim[2] - size[2]);
  s_in_offset  =  element * (dim[2] * (dim[1] - size[1]));
  m_out_offset = (element * size[0]) - element;
  li_in        =  element * (start[2]+dim[2]*(start[1]+dim[1]*start[0]));

  for(s=0 ;s<size[0]; s++) {      /* fast changing out */
    li_out = element*s;
    for(m=0; m<size[1]; m++) {    /* slow changing out */
      for(f=0 ;f<size[2]; f++) {  /* mid  changing out */
	for(e=0; e<element; e++) out[li_out++] = in[li_in++];
	li_out += m_out_offset;
      }
      li_in  += m_in_offset;
    }
    li_in += s_in_offset;
  }

}

void fft_pack_block_permute2(double *in, double *out, int start[3], int size[3], 
			 int dim[3],int element)
{
  /* slow,mid and fast changing indices for input  grid */
  int s,m,f,e;
  /* linear index of in grid, linear index of out grid */
  int li_in,li_out=0;
  /* offsets for indizes in input grid */
  int m_in_offset,s_in_offset;
  /* offset for slow changing index of output grid */
  int s_out_offset;
  /* start index for mid changing index of output grid */
  int m_out_start;

  m_in_offset  = element * (dim[2]-size[2]);
  s_in_offset  = element * (dim[2] * (dim[1]-size[1]));
  s_out_offset = (element * size[0] * size[1]) - element;
  li_in        = element * (start[2]+dim[2]*(start[1]+dim[1]*start[0]));

  for(s=0 ;s<size[0]; s++) {      /* mid changing out */
    m_out_start = element*(s * size[1]);
    for(m=0; m<size[1]; m++) {    /* fast changing out */
      li_out = m_out_start + element*m;
      for(f=0 ;f<size[2]; f++) {  /* slow  changing out */
	for(e=0; e<element; e++) out[li_out++] = in[li_in++];
	li_out += s_out_offset;
      }
      li_in += m_in_offset; 
    }
    li_in += s_in_offset; 
  }

}

void fft_unpack_block(double *in, double *out, int start[3], int size[3], 
		  int dim[3], int element)
{
  /* mid and slow changing indices */
  int m,s;
  /* linear index of in grid, linear index of out grid */
  int li_in=0,li_out;
  /* copy size */
  int copy_size;
  /* offset for indizes in input grid */
  int m_in_offset;
  /* offsets for indizes in output grid */
  int m_out_offset,s_out_offset;

  copy_size    = element * (size[2] * sizeof(double));
  m_out_offset = element * dim[2];
  s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  m_in_offset  = element * size[2];
  li_out       = element * (start[2]+dim[2]*(start[1]+dim[1]*start[0]));

  for(s=0 ;s<size[0]; s++) {
    for(m=0; m<size[1]; m++) {
      memcpy(&(out[li_out]), &(in[li_in]), copy_size);
      li_in  += m_in_offset;
      li_out += m_out_offset;
    }
    li_out += s_out_offset;
  }

}

/************************************************
 * privat functions
 ************************************************/

int fft_find_comm_groups(int grid1[3], int grid2[3], int *node_list1, int *node_list2, 
		     int *group, int *pos, int *my_pos)
{
  int i;
  /* communication group cell size on grid1 and grid2 */
  int s1[3], s2[3];
  /* The communication group cells build the same super grid on grid1 and grid2 */
  int ds[3];
  /* communication group size */
  int g_size=1;
  /* comm. group cell index */
  int gi[3];
  /* position of a node in a grid */
  int p1[3], p2[3];
  /* node identity */
  int n;
  /* this_node position in the communication group. */
  int c_pos=-1;
  /* flag for group identification */
  int my_group=0;

  FFT_TRACE(fprintf(stderr,"%d: fft_find_comm_groups:\n",this_node));
  FFT_TRACE(fprintf(stderr,"%d: for grid1=(%d,%d,%d) and grids=(%d,%d,%d)\n",
		    this_node,grid1[0],grid1[1],grid1[2],grid2[0],grid2[1],grid2[2]));

  /* calculate dimension of comm. group cells for both grids */ 
  if( (grid1[0]*grid1[1]*grid1[2]) != (grid2[0]*grid2[1]*grid2[2]) ) return -1; /* unlike number of nodes */
  for(i=0;i<3;i++) {
    s1[i] = grid1[i] / grid2[i];
    if(s1[i] == 0) s1[i] = 1;
    else if(grid1[i] != grid2[i]*s1[i]) return -1; /* grids do not match!!! */

    s2[i] = grid2[i] / grid1[i];
    if(s2[i] == 0) s2[i] = 1;
    else if(grid2[i] != grid1[i]*s2[i]) return -1; /* grids do not match!!! */

    ds[i] = grid2[i] / s2[i]; 
    g_size *= s2[i];
  }

  /* calc node_list2 */
  /* loop through all comm. group cells */
  for(gi[2] = 0; gi[2] < ds[2]; gi[2]++) 
    for(gi[1] = 0; gi[1] < ds[1]; gi[1]++)
      for(gi[0] = 0; gi[0] < ds[0]; gi[0]++) {
	/* loop through all nodes in that comm. group cell */
	for(i=0;i<g_size;i++) {
	  p1[0] = (gi[0]*s1[0]) + (i%s1[0]);
	  p1[1] = (gi[1]*s1[1]) + ((i/s1[0])%s1[1]);
	  p1[2] = (gi[2]*s1[2]) + (i/(s1[0]*s1[1]));

	  p2[0] = (gi[0]*s2[0]) + (i%s2[0]);
	  p2[1] = (gi[1]*s2[1]) + ((i/s2[0])%s2[1]);
	  p2[2] = (gi[2]*s2[2]) + (i/(s2[0]*s2[1]));

	  n = node_list1[ get_linear_index(p1[0],p1[1],p1[2],grid1) ];
	  node_list2[ get_linear_index(p2[0],p2[1],p2[2],grid2) ] = n ;

	  pos[3*n+0] = p2[0];  pos[3*n+1] = p2[1];  pos[3*n+2] = p2[2];	  
	  if(my_group==1) group[i] = n;
	  if(n==this_node && my_group==0) { 
	    my_group = 1; 
	    c_pos = i;
	    my_pos[0] = p2[0]; my_pos[1] = p2[1]; my_pos[2] = p2[2];
	    i=-1; /* restart the loop */ 
	  }
	}
	my_group=0;
      }

  /* permute comm. group according to the nodes position in the group */
  /* This is necessary to have matching node pairs during communication! */
  while( c_pos>0 ) {
    n=group[g_size-1];
    for(i=g_size-1; i>0; i--) group[i] = group[i-1];
    group[0] = n;
    c_pos--;
  }
  return g_size;
}

int fft_calc_local_mesh(int n_pos[3], int n_grid[3], int mesh[3], double mesh_off[3], 
		     int loc_mesh[3], int start[3])
{
  int i, last[3], size=1;
  
  for(i=0;i<3;i++) {
    start[i] = (int)ceil( (mesh[i]/(double)n_grid[i])*n_pos[i]     - mesh_off[i] );
    last[i]  = (int)floor((mesh[i]/(double)n_grid[i])*(n_pos[i]+1) - mesh_off[i] );
    /* correct round off errors */
    if( (mesh[i]/(double)n_grid[i])*(n_pos[i]+1) - mesh_off[i] - last[i] < 1.0e-15 ) last[i]--;
    if(1.0+ (mesh[i]/(double)n_grid[i])*n_pos[i]-mesh_off[i]-start[i] < 1.0e-15 ) start[i]--;
    loc_mesh[i] = last[i]-start[i]+1;
    size *= loc_mesh[i];
  }
  return size;
}


int fft_calc_send_block(int pos1[3], int grid1[3], int pos2[3], int grid2[3], 
		    int mesh[3], double mesh_off[3], int block[6])
{
  int i,size=1;
  int mesh1[3], first1[3], last1[3];
  int mesh2[3], first2[3], last2[3];

  fft_calc_local_mesh(pos1, grid1, mesh, mesh_off, mesh1, first1);
  fft_calc_local_mesh(pos2, grid2, mesh, mesh_off, mesh2, first2);

  for(i=0;i<3;i++) {
    last1[i] = first1[i] + mesh1[i] -1;
    last2[i] = first2[i] + mesh2[i] -1;
    block[i  ] = imax(first1[i],first2[i]) - first1[i];
    block[i+3] = (imin(last1[i], last2[i] ) - first1[i])-block[i]+1;
    size *= block[i+3];
  }
  return size;
}

void fft_print_fft_plan(fft_forw_plan pl)
{
  int i;

  fprintf(stderr,"%d: dir=%d, row_dir=%d, n_permute=%d, n_ffts=%d\n",
	  this_node, pl.dir,  pl.row_dir, pl.n_permute, pl.n_ffts);

  fprintf(stderr,"%d:    local: old_mesh=(%d,%d,%d), new_mesh=(%d,%d,%d), start=(%d,%d,%d)\n",this_node,
	  pl.old_mesh[0],  pl.old_mesh[1],  pl.old_mesh[2], 
	  pl.new_mesh[0],  pl.new_mesh[1],  pl.new_mesh[2], 
	  pl.start[0], pl.start[1],  pl.start[2]);

  fprintf(stderr,"%d:    g_size=%d group=(",this_node,pl.g_size);
  for(i=0;i<pl.g_size-1;i++) fprintf(stderr,"%d,", pl.group[i]);
  fprintf(stderr,"%d)\n",pl.group[pl.g_size-1]);

  fprintf(stderr,"%d:    send=[",this_node);
  for(i=0;i<pl.g_size;i++) fprintf(stderr,"(%d,%d,%d)+(%d,%d,%d), ",
				   pl.send_block[6*i+0], pl.send_block[6*i+1], pl.send_block[6*i+2],
				   pl.send_block[6*i+3], pl.send_block[6*i+4], pl.send_block[6*i+5]);
  fprintf(stderr,"]\n%d:    recv=[",this_node);
  for(i=0;i<pl.g_size;i++) fprintf(stderr,"(%d,%d,%d)+(%d,%d,%d), ",
				   pl.recv_block[6*i+0], pl.recv_block[6*i+1], pl.recv_block[6*i+2],
				   pl.recv_block[6*i+3], pl.recv_block[6*i+4], pl.recv_block[6*i+5]);
  fprintf(stderr,"]\n");
 
 fflush(stderr);
}

void fft_print_global_fft_mesh(fft_forw_plan plan, double *data, int element, int num)
{
  int i0,i1,i2,b=1;
  int mesh,divide=0,block1=-1,start1;
  int st[3],en[3],si[3];
  int my=-1;
  double tmp;

  for(i1=0;i1<3;i1++) {
    st[i1] = plan.start[i1];
    en[i1] = plan.start[i1]+plan.new_mesh[i1];
    si[i1] = plan.new_mesh[i1];
  }

  mesh = plan.new_mesh[2];
  MPI_Barrier(comm_cart);  
  if(this_node==0) fprintf(stderr,"All: Print Global Mesh: (%d of %d elements)\n",
			   num+1,element);
  MPI_Barrier(comm_cart);
  for(i0=0;i0<n_nodes;i0++) {
    MPI_Barrier(comm_cart);
    if(i0==this_node) fprintf(stderr,"%d: range (%d,%d,%d)-(%d,%d,%d)\n",this_node,st[0],st[1],st[2],en[0],en[1],en[2]);
  }
  MPI_Barrier(comm_cart);
  while(divide==0) {
    if(b*mesh > 7) {
      block1=b;
      divide = (int)ceil(mesh/(double)block1);
    }
    b++;
  }

  for(b=0;b<divide;b++) {
    start1 = b*block1;
    for(i0=mesh-1; i0>=0; i0--) {
      for(i1=start1; i1<imin(start1+block1,mesh);i1++) {
	for(i2=0; i2<mesh;i2++) {
	  if(i0>=st[0] && i0<en[0] && i1>=st[1] && 
	     i1<en[1] && i2>=st[2] && i2<en[2]) my=1;
	  else my=0;
	  MPI_Barrier(comm_cart);
	  if(my==1) {
	   
	    tmp=data[num+(element*((i2-st[2])+si[2]*((i1-st[1])+si[1]*(i0-st[0]))))];
	    if(fabs(tmp)>1.0e-15) {
	      if(tmp<0) fprintf(stderr,"%1.2e",tmp);
	      else      fprintf(stderr," %1.2e",tmp);
	    }
	    else {
	      fprintf(stderr," %1.2e",0.0);
	    }
	  }
	  MPI_Barrier(comm_cart);
	}
	if(my==1) fprintf(stderr," | ");
      }
      if(my==1) fprintf(stderr,"\n");
    }
    if(my==1) fprintf(stderr,"\n");
  }

}


#endif /* defined(P3M) || defined(DP3M) */
