#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*#include <dxmldef.h> */
#include "global.h"
#include "cells.h"
#include "grid.h"
#include "particle_data.h"
#include "utils.h"
#include "global.h"
#include "particle_data.h"
#include "myp3m.h"
#include "communication.h"

p3m_struct p3m;
int **meshift;

/* stores the interpolation of the charge assignment functions */
int MaxInterpol = 4*50048;
double  **intCAF; 
/* charge of the charged particles */
double  *QP;

/* mesh contributions of the charged particles */
double  *QL;
/* global index of a charged particle */
int     *global;
/* ...and for the FFT: */
double  ***Q_re;
double  ***Q_im;
double  ***Q_re_tot;
double modadd1, modadd2;

/* tmp variables*/
int local_mesh[3];

double Hi;
double inner_left_gridpoint[3], inner_right_gridpoint[3];/* coordinates of local grid*/
double left_gridpoint[3], right_gridpoint[3];/*coordinates of ext grid*/
double ca_margin[6];/* margin of gridpoints to local box, in units of 1/Hi */
int send_number[6], rec_number[6];
int grid_margin[6];/* margin of gridpoints around local mesh*/
double ld[6][3], ur[6][3];/* points left down & right up for send boxes*/
double rld[6][3], rur[6][3];/* same for recieve boxes*/
int ext_mesh[3];/* extended mesh, ie local mesh + grid_margin*/
double *send_grids[6];/* pointer tb send*/

/*
int zfft_apply_3d_( char *, char *, char *, ... );
int zfft_init_3d_( int *, int *, int *, DXML_Z_FFT_STRUCTURE_3D *, int * );
int zfft_exit_3d_( DXML_Z_FFT_STRUCTURE_3D * ); 
old dxml routines
*/

void setp3m();   
void   P3M_init();
void   compute_measures();
void   P3M_perform();
void   P3M_exit();
void interpolate_charge_assignment_function(void);
void calculate_meshift(void);


void setp3m()
{
  /* just for the building process, should happen in some script*/
p3m.P = 7;
p3m.alpha = 0.2;
p3m.mesh[0] = 150;
p3m.mesh[1] = 150;
p3m.mesh[2] = 150;
printf("%d setp3m finished\n",this_node);

}
/*----------------------------------------------------------------------*/
void   P3M_init()
{
  int i,j,k;
  setp3m();
  /* compute local mesh for fixed processor box */
  /* Direction of grid indicees in coordinate 2 is inconsistent!!!!!*/
  Hi = (double)(p3m.mesh[0]-1)/box_l[0];

  compute_measures();
  intCAF = malloc(sizeof(double *)*(p3m.P));
  for(i=0;i<p3m.P;i++)
    intCAF[i] = malloc(sizeof(double)*(2*MaxInterpol+1));

  meshift = malloc(sizeof(double*)*3);
  for(i=0;i<3;i++)
    meshift[i] = malloc(sizeof(double)*p3m.mesh[i]);
  Q_re = malloc(sizeof(double *)*ext_mesh[0]);
  Q_re_tot = malloc(sizeof(double *)*ext_mesh[0]);
  Q_im = malloc(sizeof(double *)*ext_mesh[0]);
  for(i=0;i<ext_mesh[0];i++)
    {
      Q_re[i] = malloc(sizeof(double *)*ext_mesh[1]);
      Q_re_tot[i] = malloc(sizeof(double *)*ext_mesh[1]);
      Q_im[i] = malloc(sizeof(double *)*ext_mesh[1]);
      for(j=0;j<ext_mesh[1];j++)
	{
	  Q_re[i][j] = malloc(sizeof(double)*ext_mesh[2]);
	  Q_re_tot[i][j] = malloc(sizeof(double)*ext_mesh[2]);
	  Q_im[i][j] = malloc(sizeof(double)*ext_mesh[2]);
	}
    }

  interpolate_charge_assignment_function();
  calculate_meshift();

  switch (p3m.P) {
  case 2 : case 4 : case 6 : 
    { modadd1 = 0.5; modadd2 = -0.5;} break;
  case 1 :case 3 : case 5 : case 7 : 
    { modadd1 = 0.0; modadd2 =  0.5;} break;
  default : {
    fprintf(stderr,"Error in function 'P3M_init':\n");
    fprintf(stderr,"Charge assignment order P=%d unknown.\nProgram terminated\n\n",p3m.P);
    exit(1);
  } break;
  }
printf("%d P3M_init finished\n",this_node);

}
/*----------------------------------------------------------------------*/
void compute_measures()
{
  int i,j,k;
  for(i=0;i<3;i++)
    {
      local_mesh[i] = (int) ((local_box_l[i]+2.*skin)*Hi);
      inner_left_gridpoint[i]  = ceil( (my_left[i]-skin) *Hi)/Hi - (my_left[i]-skin);
      inner_right_gridpoint[i] = floor((my_right[i]+skin)*Hi)/Hi - my_right[i]+skin;
    }
  printf("loc: %f %f %f\n",local_box_l[0],skin,Hi);
  printf("%d\n",local_mesh[0]);
  for(i=0;i<6;i++)
    {
      if(i<=1) j=0;
      else if(i<=3) j=2;
      else j=1;
      if(i%2 == 0) ca_margin[i]   =  ( ceil((my_left[i] -skin)*Hi)/Hi - (my_left[i]-skin))/Hi;
      else         ca_margin[i]   =  (floor((my_right[i]+skin)*Hi)/Hi -  my_right[i]+skin)/Hi;
      /*ca_margin[i] += skin/Hi;*/
    }
  for(i=0;i<3;i++)
    {
      if(i==0) {j=0;k=1;}
      else if(i==1) {j=4;k=5;}
      else {j=2;k=3;}
      left_gridpoint[i]  = inner_left_gridpoint[i]  - (((double) p3m.P/2. + ca_margin[j])  )/Hi;
      grid_margin[j] = (int) ((double) p3m.P/2. + ca_margin[j]);
      right_gridpoint[i] = inner_right_gridpoint[i] + (((double) p3m.P/2. + ca_margin[k]))/Hi;
      grid_margin[k] = (int) ((double) p3m.P/2. + ca_margin[k]);
      ext_mesh[i] = local_mesh[i] + grid_margin[j] + grid_margin[k];
    }
  printf("thisnode: %d\n",this_node);
  printf("local_grid: %d %d %d \n",local_mesh[0],local_mesh[1], local_mesh[2]);
  printf("gridmargin: %d %d %d %d %d %d \n",grid_margin[0],grid_margin[1],grid_margin[2],
	 grid_margin[3],grid_margin[4],grid_margin[5]);
  printf("ca_margin: %f\t%f\t%f\n%f\t%f\t%f\n",ca_margin[0],ca_margin[1],ca_margin[2],
	 ca_margin[3],ca_margin[4],ca_margin[5]);

  exit(0);

  ld[0][0]       =  0;
  ld[0][1]       =  0;
  ld[0][2]       =  ext_mesh[2];
  ur[0][0]       =  grid_margin[0];
  ur[0][1]       =  ext_mesh[1];
  ur[0][2]       =  0;

  ld[1][0]       =  grid_margin[0]+local_mesh[0];
  ld[1][1]       =  0;
  ld[1][2]       =  ext_mesh[2];
  ur[1][0]       =  ext_mesh[0];
  ur[1][1]       =  ext_mesh[1];
  ur[1][2]       =  0;

  ld[2][0]       =  grid_margin[0];
  ld[2][1]       =  0;
  ld[2][2]       =  ext_mesh[2];
  ur[2][0]       =  grid_margin[0]+local_mesh[0];
  ur[2][1]       =  ext_mesh[1];
  ur[2][2]       =  grid_margin[3]+local_mesh[2];

  ld[3][0]       =  grid_margin[0];
  ld[3][1]       =  0;
  ld[3][2]       =  grid_margin[3];
  ur[3][0]       =  grid_margin[0]+local_mesh[0];
  ur[3][1]       =  ext_mesh[1];
  ur[3][2]       =  0;

  ld[4][0]       =  grid_margin[0];
  ld[4][1]       =  0;
  ld[4][2]       =  grid_margin[3]+local_mesh[2];
  ur[4][0]       =  grid_margin[0]+local_mesh[0];
  ur[4][1]       =  grid_margin[4];
  ur[4][2]       =  grid_margin[3];

  ld[5][0]       =  grid_margin[0];
  ld[5][1]       =  grid_margin[4]+local_mesh[1];
  ld[5][2]       =  grid_margin[3]+local_mesh[2];
  ur[5][0]       =  grid_margin[0]+local_mesh[0];
  ur[5][1]       =  ext_mesh[1];
  ur[5][2]       =  grid_margin[3];

  send_number[0] =  ext_mesh[1]   * ext_mesh[2]   * grid_margin[0];
  send_number[1] =  ext_mesh[1]   * ext_mesh[2]   * grid_margin[1];
  send_number[2] =  local_mesh[0] * ext_mesh[1]   * grid_margin[2];
  send_number[3] =  local_mesh[0] * ext_mesh[1]   * grid_margin[3];
  send_number[4] =  local_mesh[0] * local_mesh[2] * grid_margin[4];
  send_number[5] =  local_mesh[0] * local_mesh[2] * grid_margin[5];

  for(i=0;i<6;i++)
    send_grids[i] = malloc(sizeof(double)*send_number[i]);

  rec_number[0] = (local_mesh[1]+(int)((double)p3m.P/2+1.-ca_margin[4]) +(int)((double)p3m.P/2+1.-ca_margin[5]))*
                  (local_mesh[2]+(int)((double)p3m.P/2+1.-ca_margin[2]) +(int)((double)p3m.P/2+1.-ca_margin[3]))*
                  (int)((double)p3m.P/2+1.-ca_margin[0]);

  rec_number[1] = (local_mesh[1]+(int)((double)p3m.P/2+1.-ca_margin[4]) +(int)((double)p3m.P/2+1.-ca_margin[5]))*
                  (local_mesh[2]+(int)((double)p3m.P/2+1.-ca_margin[2]) +(int)((double)p3m.P/2+1.-ca_margin[3]))*
                  (int)((double)p3m.P/2+1.-ca_margin[1]);

  rec_number[2] =  local_mesh[0]*
                  (local_mesh[1]+(int)((double)p3m.P/2+1.-ca_margin[4]) +(int)((double)p3m.P/2+1.-ca_margin[5]))*
                  (int)((double)p3m.P/2+1.-ca_margin[2]);

  rec_number[3] =  local_mesh[0]*
                  (local_mesh[1]+(int)((double)p3m.P/2+1.-ca_margin[4]) +(int)((double)p3m.P/2+1.-ca_margin[5]))*
                  (int)((double)p3m.P/2+1.-ca_margin[3]);

  rec_number[4] = local_mesh[0] * local_mesh[2] * (int)((double)p3m.P/2+1.-ca_margin[4]);

  rec_number[5] = local_mesh[0] * local_mesh[2] * (int)((double)p3m.P/2+1.-ca_margin[5]);

  rld[0][0] = ld[0][0] + grid_margin[0];
  rld[0][1] = ld[0][1];
  rld[0][2] = ld[0][2];
  rur[0][0] = ur[0][0] + grid_margin[0];
  rur[0][1] = ur[0][1];
  rur[0][2] = ur[0][2];

  rld[1][0] = ld[1][0] - grid_margin[1];
  rld[1][1] = ld[1][1];
  rld[1][2] = ld[1][2];
  rur[1][0] = ur[1][0] - grid_margin[1];
  rur[1][1] = ur[1][1];
  rur[1][2] = ur[1][2];


  rld[2][0] = ld[2][0];
  rld[2][1] = ld[2][1];
  rld[2][2] = ld[2][2] - grid_margin[2];
  rur[2][0] = ur[2][0];
  rur[2][1] = ur[2][1];
  rur[2][2] = ur[2][2] - grid_margin[2];


  rld[3][0] = ld[3][0];
  rld[3][1] = ld[3][1];
  rld[3][2] = ld[3][2] + grid_margin[3];
  rur[3][0] = ur[3][0];
  rur[3][1] = ur[3][1];
  rur[3][2] = ur[3][2] + grid_margin[3];


  rld[4][0] = ld[4][0];
  rld[4][1] = ld[4][1] + grid_margin[4];
  rld[4][2] = ld[4][2];
  rur[4][0] = ur[4][0];
  rur[4][1] = ur[4][1] + grid_margin[4];
  rur[4][2] = ur[4][2];

  rld[5][0] = ld[5][0];
  rld[5][1] = ld[5][1] - grid_margin[5];
  rld[5][2] = ld[5][2];
  rur[5][0] = ur[5][0];
  rur[5][1] = ur[5][1] - grid_margin[5];
  rur[5][2] = ur[5][2];


}
/*----------------------------------------------------------------------*/
void P3M_perform()
{
  int i,j,k,l,m,n;
  int QZahl = 0;
  int assignshift[3];
  

  double d1;
  double **G;
  int xarg, yarg, zarg;
  int xpos, ypos, zpos;

  double MI2 = 2.0*(double)MaxInterpol;
  int Gi0, Gi1, Gi2;
  double T1, T2, T3;



  assignshift[0] = local_mesh[0]-(p3m.P-1)/2;  /* INTEGER division is vital here! */
  assignshift[1] = local_mesh[1]-(p3m.P-1)/2; 
  assignshift[2] = local_mesh[2]-(p3m.P-1)/2;
  QP = malloc(sizeof(double)*n_particles);
  QL = malloc(sizeof(double)*n_particles*(p3m.P*p3m.P*p3m.P));
  global = malloc(sizeof(int)*n_particles);
  for (i=0; i<n_particles; i++) if (fabs(particles[i].q)>1e-5) 
    {
      global[QZahl] = i;
      QP[QZahl++] = particles[i].q;
    }
  G = malloc(sizeof(double*)*QZahl);
  for(i=0;i<QZahl;i++)
    G[i] = malloc(sizeof(double *)*3);

  /* init Q_re and Q_im */
  for (i=0; i<ext_mesh[0]; i++) 
    for (j=0; j<ext_mesh[1]; j++) 
      for (k=0; k<ext_mesh[2]; k++)
	{
	  Q_re[i][j][k] = Q_im[i][j][k] = 0.0;
	  Q_re_tot[i][j][k] = Q_im[i][j][k] = 0.0;
	}

  /*  m = 0;*/ /* <-- enumerate the charge fractions consecutively */
  /* G: first grid for the charge smearing */
  for (i=0; i<QZahl; i++) 
    {
      d1      = (particles[global[i]].p[0]-my_left[0])*Hi + modadd1; 
      G[i][0] = Gi0 = (int)(d1 + modadd2) + assignshift[0];
      xarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      d1      = (particles[global[i]].p[1]-my_left[1])*Hi + modadd1; 
      G[i][1] = Gi1 = (int)(d1 + modadd2) + assignshift[1];
      yarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );

      d1      = (particles[global[i]].p[2]-my_left[2])*Hi + modadd1; 
      G[i][2] = Gi2 = (int)(d1 + modadd2) + assignshift[2];
      zarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );

      for (j = 0; j < p3m.P; j++) {
	xpos = (Gi0 + j) % local_mesh[0];
	T1 = QP[i] * intCAF[j][xarg];
	for (k = 0; k < p3m.P; k++) {
	  ypos = (Gi1 + k) % local_mesh[1];
	  T2 = T1 * intCAF[k][yarg];
	  for (l = 0; l < p3m.P; l++) {
	    zpos = (Gi2 + l) % local_mesh[2];
	    T3 = T2 * intCAF[l][zarg];
	if(xpos<0) exit(1);
	    Q_re[xpos][ypos][zpos] += T3;
	  }
	}
      }
    }

  m=0;
  for(n=0;n<6;n++)
    {
      /* here: update grids*/
      for(i=ld[n][0];i<ur[n][0];i++)
	for(j=ld[n][1];j<ur[n][1];j++)
	  for(k=ld[n][2];k<ur[n][2];k++)
	    send_grids[n][m++] = Q_re[i][j][k];
      /*here: communication*/
    }
    
}
/*----------------------------------------------------------------------*/
void   P3M_exit()
{
  free(intCAF);
  free(meshift);
}
/*----------------------------------------------------------------------*/
void interpolate_charge_assignment_function(void)
{
  /* Interpolates the P-th order charge assignment function from
     Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
     are also tabulated in Deserno/Holm. */

  double dInterpol=(double)MaxInterpol, x;
  long   i;

  fprintf(stderr,"%d - interpolating the order-%d charge assignment function",this_node,p3m.P);

  switch (p3m.P) {
  case 1 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 1.0;
    }
  } break;
  case 2 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 0.5-x;
      intCAF[1][i+MaxInterpol] = 0.5+x;
    }
  } break;
  case 3 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 0.5*SQR(0.5 - x);
      intCAF[1][i+MaxInterpol] = 0.75 - SQR(x);
      intCAF[2][i+MaxInterpol] = 0.5*SQR(0.5 + x);
    }
  } break;
  case 4 :{ 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
      intCAF[1][i+MaxInterpol] = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
      intCAF[2][i+MaxInterpol] = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
      intCAF[3][i+MaxInterpol] = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    }
  } break;
  case 5 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
      intCAF[1][i+MaxInterpol] = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
      intCAF[2][i+MaxInterpol] = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
      intCAF[3][i+MaxInterpol] = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
      intCAF[4][i+MaxInterpol] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    }
  } break;
  case 6 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
      intCAF[1][i+MaxInterpol] = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
      intCAF[2][i+MaxInterpol] = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
      intCAF[3][i+MaxInterpol] = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
      intCAF[4][i+MaxInterpol] = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
      intCAF[5][i+MaxInterpol] = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    }
  } break;
  case 7 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
      intCAF[1][i+MaxInterpol] = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
      intCAF[2][i+MaxInterpol] = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
      intCAF[3][i+MaxInterpol] = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
      intCAF[4][i+MaxInterpol] = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
      intCAF[5][i+MaxInterpol] = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
      intCAF[6][i+MaxInterpol] = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    } 
  } break;
  default :{
    fprintf(stderr,"Error in function 'interpolate_charge_assignment_function':");
    fprintf(stderr,"Charge assignment order %d unknown.\nProgram terminated.\n\n",p3m.P);
    exit(1);
  }
  }

  fprintf(stderr,"\n");
}

/*----------------------------------------------------------------------*/
void calculate_meshift(void)
{
  /* shifts the mesh points by mesh/2 */
  
  int i,j;

  fprintf(stderr," - calculating mesh-shift");
  
  for (i=0; i<3; i++)
    for (j=0; j<p3m.mesh[i]; j++) 
      meshift[i][j] = j - dround(j/ (double) p3m.mesh[i])*(double) p3m.mesh[i]; 
  
  fprintf(stderr,"\n");
}
/*----------------------------------------------------------------------*/
