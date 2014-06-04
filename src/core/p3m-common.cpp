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
/** \file p3m-common.cpp P3M main file.
*/
#include "p3m-common.hpp"

#if defined(P3M) || defined(DP3M)

void p3m_common_parameter_pre_init(p3m_parameter_struct *params) {
  params->alpha = 0.0;
  params->alpha_L = 0.0;
  params->r_cut = 0.0;
  params->r_cut_iL = 0.0;
  params->mesh[0] = 0;
  params->mesh[1] = 0;
  params->mesh[2] = 0;
  params->cao = 0;
  params->cao3 = 0;
  params->mesh_off[0] = P3M_MESHOFF;
  params->mesh_off[1] = P3M_MESHOFF;
  params->mesh_off[2] = P3M_MESHOFF;
  params->inter = P3M_N_INTERPOL;
  params->inter2 = 0;
  params->accuracy = 0.0;
  params->epsilon = P3M_EPSILON;
  params->cao_cut[0] = 0.0;
  params->cao_cut[1] = 0.0;
  params->cao_cut[2] = 0.0;
  params->a[0] = 0.0;
  params->a[1] = 0.0;
  params->a[2] = 0.0;
  params->ai[0] = 0.0;
  params->ai[1] = 0.0;
  params->ai[2] = 0.0;
  params->additional_mesh[0] = 0;
  params->additional_mesh[1] = 0;
  params->additional_mesh[2] = 0;
}

/** Debug function printing p3m structures */
void p3m_p3m_print_local_mesh(p3m_local_mesh l) 
{
  fprintf(stderr,"%d: p3m_local_mesh: dim=(%d,%d,%d), size=%d\n",this_node,
	  l.dim[0],l.dim[1],l.dim[2],l.size);
  fprintf(stderr,"%d:    ld_ind=(%d,%d,%d), ld_pos=(%f,%f,%f)\n",this_node,
	  l.ld_ind[0],l.ld_ind[1],l.ld_ind[2],
	  l.ld_pos[0],l.ld_pos[1],l.ld_pos[2]);
  fprintf(stderr,"%d:    inner=(%d,%d,%d) [(%d,%d,%d)-(%d,%d,%d)]\n",this_node,
	  l.inner[0],l.inner[1],l.inner[2],
	  l.in_ld[0],l.in_ld[1],l.in_ld[2],
	  l.in_ur[0],l.in_ur[1],l.in_ur[2]);
  fprintf(stderr,"%d:    margin = (%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.margin[0],l.margin[1],l.margin[2],l.margin[3],l.margin[4],l.margin[5]);
  fprintf(stderr,"%d:    r_margin=(%d,%d, %d,%d, %d,%d)\n",this_node,
	  l.r_margin[0],l.r_margin[1],l.r_margin[2],l.r_margin[3],l.r_margin[4],l.r_margin[5]);
}

/** Debug function printing p3m structures */
void p3m_p3m_print_send_mesh(p3m_send_mesh sm) 
{
  int i;
  fprintf(stderr,"%d: p3m_send_mesh: max=%d\n",this_node,sm.max);
  for(i=0;i<6;i++) {
    fprintf(stderr,"%d:  dir=%d: s_dim (%d,%d,%d)  s_ld (%d,%d,%d) s_ur (%d,%d,%d) s_size=%d\n",this_node,i,sm.s_dim[i][0],sm.s_dim[i][1],sm.s_dim[i][2],sm.s_ld[i][0],sm.s_ld[i][1],sm.s_ld[i][2],sm.s_ur[i][0],sm.s_ur[i][1],sm.s_ur[i][2],sm.s_size[i]);
    fprintf(stderr,"%d:         r_dim (%d,%d,%d)  r_ld (%d,%d,%d) r_ur (%d,%d,%d) r_size=%d\n",this_node,sm.r_dim[i][0],sm.r_dim[i][1],sm.r_dim[i][2],sm.r_ld[i][0],sm.r_ld[i][1],sm.r_ld[i][2],sm.r_ur[i][0],sm.r_ur[i][1],sm.r_ur[i][2],sm.r_size[i]);
  }
}

void p3m_add_block(double *in, double *out, int start[3], int size[3], int dim[3])
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

double p3m_analytic_cotangent_sum(int n, double mesh_i, int cao)
{
  double c, res=0.0;
  c = SQR(cos(PI*mesh_i*(double)n));

  switch (cao) {
  case 1 : { 
    res = 1; 
    break; }
  case 2 : { 
    res = (1.0+c*2.0)/3.0; 
    break; }
  case 3 : { 
    res = (2.0+c*(11.0+c*2.0))/15.0; 
    break; }
  case 4 : { 
    res = (17.0+c*(180.0+c*(114.0+c*4.0)))/315.0; 
    break; }
  case 5 : { 
    res = (62.0+c*(1072.0+c*(1452.0+c*(247.0+c*2.0))))/2835.0; 
    break; }
  case 6 : { 
    res = (1382.0+c*(35396.0+c*(83021.0+c*(34096.0+c*(2026.0+c*4.0)))))/155925.0; 
    break; }
  case 7 : { 
    res = (21844.0+c*(776661.0+c*(2801040.0+c*(2123860.0+c*(349500.0+c*(8166.0+c*4.0))))))/6081075.0; 
    break; }
  default : {
    fprintf(stderr,"%d: INTERNAL_ERROR: The value %d for the interpolation order should not occur!\n",this_node, cao);
    errexit();
  }
  }
  
  return res;
}

/** Computes the  assignment function of for the \a i'th degree
    at value \a x. */
double p3m_caf(int i, double x, int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 5 : {
    switch (i) {
    case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
    case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
    case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
    case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
    case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 6 : {
    switch (i) {
    case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
    case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
    case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
    case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
    case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
    case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  case 7 : {
    switch (i) {
    case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
    case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
    case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
    case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
    case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
    case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
    case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    default:
      fprintf(stderr,"%d: Tried to access charge assignment function of degree %d in scheme of order %d.\n",this_node,i,cao_value);
      return 0.0;
    }
  }
  default :{
    fprintf(stderr,"%d: Charge assignment order %d unknown.\n",this_node,cao_value);
    return 0.0;
  }}}
}
#endif /* defined(P3M) || defined(DP3M) */
