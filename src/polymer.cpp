/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file polymer.cpp
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and crosslinking the polymers if necessary.
 
    The corresponding header file is polymer.hpp.
 
    Created:       27.02.2003 by BAM
       Based upon 'polymer.tcl' by BAM (20.02.2003).
*/

#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <cmath>
#include "utils.hpp"
#include "polymer.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "interaction_data.hpp"
#include "random.hpp"
#include "integrate.hpp"
#include "constraint.hpp"




/************************************************************* 
 * Functions                                                 *
 * ---------                                                 *
 *************************************************************/



int mindist3(int part_id, double r_catch, int *ids) {
  Particle *partCfgMD;
  double dx,dy,dz;
  int i, me, caught=0;

  partCfgMD = (Particle*)malloc(n_part*sizeof(Particle));
  mpi_get_particles(partCfgMD, NULL);
  me = -1; /* Since 'mpi_get_particles' returns the particles unsorted, it's most likely that 'partCfgMD[i].p.identity != i'
	      --> prevent that! */
  for(i=0; i<n_part; i++) if (partCfgMD[i].p.identity == part_id) me = i; 
  if (me == -1) {
    char *errtxt = runtime_error(128 + ES_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt, "{049 failed to find desired particle %d} ",part_id);
    return 0;
  }
  for (i=0; i<n_part; i++) {
    if (i != me) {
      dx = partCfgMD[me].r.p[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
      dy = partCfgMD[me].r.p[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
      dz = partCfgMD[me].r.p[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
      if (sqrt(SQR(dx)+SQR(dy)+SQR(dz)) < r_catch) ids[caught++]=partCfgMD[i].p.identity;
    }
  }
  free(partCfgMD); 
  return (caught);
}



double mindist4(double pos[3]) {
  Particle *partCfgMD;
  double mindist=30000.0, dx,dy,dz;
  int i;

  if (n_part ==0) return (dmin(dmin(box_l[0],box_l[1]),box_l[2]));
  partCfgMD = (Particle*)malloc(n_part*sizeof(Particle));
  mpi_get_particles(partCfgMD, NULL); 
  for (i=0; i<n_part; i++) {
    dx = pos[0] - partCfgMD[i].r.p[0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = pos[1] - partCfgMD[i].r.p[1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = pos[2] - partCfgMD[i].r.p[2];   dz -= dround(dz/box_l[2])*box_l[2];
    mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
  }
  free(partCfgMD); 
  if (mindist<30000.0)
    return (sqrt(mindist));
  return (-1.0);
}

double buf_mindist4(double pos[3], int n_add, double *add) {
  double mindist=30000.0, dx,dy,dz;
  int i;

  if (n_add == 0) return (dmin(dmin(box_l[0],box_l[1]),box_l[2]));
  for (i=0; i<n_add; i++) {
    dx = pos[0] - add[3*i + 0];   dx -= dround(dx/box_l[0])*box_l[0];
    dy = pos[1] - add[3*i + 1];   dy -= dround(dy/box_l[1])*box_l[1];
    dz = pos[2] - add[3*i + 2];   dz -= dround(dz/box_l[2])*box_l[2];
    mindist = dmin(mindist, SQR(dx)+SQR(dy)+SQR(dz));
  }
  if (mindist<30000.0) return (sqrt(mindist));
  return (-1.0);
}



int collision(double pos[3], double shield, int n_add, double *add) {
  if (mindist4(pos) > shield && buf_mindist4(pos, n_add, add) > shield) return (0);
  return (1);
}

#ifdef CONSTRAINTS

int constraint_collision(double *p1, double *p2){
  Particle part1,part2;
  double d1,d2,v[3];
  Constraint *c;
  int i;
  double folded_pos1[3];
  double folded_pos2[3];
  int img[3];

  memcpy(folded_pos1, p1, 3*sizeof(double));
  fold_position(folded_pos1, img);

  memcpy(folded_pos2, p2, 3*sizeof(double));
  fold_position(folded_pos2, img);

  for(i=0;i<n_constraints;i++){
    c=&constraints[i];
    switch(c->type){
    case CONSTRAINT_WAL:
      calculate_wall_dist(&part1,folded_pos1,&part1,&c->c.wal,&d1,v);
      calculate_wall_dist(&part2,folded_pos2,&part2,&c->c.wal,&d2,v);
      if(d1*d2<=0.0)
	return 1;
      break;
    case CONSTRAINT_SPH:
      calculate_sphere_dist(&part1,folded_pos1,&part1,&c->c.sph,&d1,v);
      calculate_sphere_dist(&part2,folded_pos2,&part2,&c->c.sph,&d2,v);
      if(d1*d2<0.0)
	return 1;
      break;
    case CONSTRAINT_CYL:
      calculate_cylinder_dist(&part1,folded_pos1,&part1,&c->c.cyl,&d1,v);
      calculate_cylinder_dist(&part2,folded_pos2,&part2,&c->c.cyl,&d2,v);
      if(d1*d2<0.0)
	return 1;
      break;
    case CONSTRAINT_MAZE:
    case CONSTRAINT_PORE:
    case CONSTRAINT_PLATE:
    case CONSTRAINT_RHOMBOID:
      break;
    }
  }
  return 0;
}

#endif

int polymerC(int N_P, int MPC, double bond_length, int part_id, double *posed, 
	     int mode, double shield, int max_try, double val_cM, int cM_dist, 
	     int type_nM, int type_cM, int type_bond, 
	     double angle, double angle2, double *posed2, int constr) {
  int p,n, cnt1,cnt2,max_cnt, bond_size, *bond, i;
  double phi,zz,rr;
  double *poly;
  double pos[3];
  double poz[3];
  double poy[3] = {0, 0, 0};
  double pox[3] = {0, 0, 0};
  double a[3] = {0, 0, 0};
  double b[3],c[3]={0., 0., 0.},d[3];
  double absc;
  poly = (double*)malloc(3*MPC*sizeof(double));

  bond_size = bonded_ia_params[type_bond].num;
  bond = (int*)malloc(sizeof(int) * (bond_size + 1));
  bond[0] = type_bond;

  cnt1 = cnt2 = max_cnt = 0;
  for (p=0; p < N_P; p++) {
    for (cnt2=0; cnt2 < max_try; cnt2++) {
      /* place start monomer */
      if (posed!=NULL) {
	/* if position of 1st monomer is given */
	if (p > 0) {
	  free(posed);
	  posed=NULL;
	} else {
	  pos[0]=posed[0];
	  pos[1]=posed[1];
	  pos[2]=posed[2];
	}
      } else {
	/* randomly set position */
	for (cnt1=0; cnt1<max_try; cnt1++) {
	  pos[0]=box_l[0]*d_random();
	  pos[1]=box_l[1]*d_random();
	  pos[2]=box_l[2]*d_random();
	  if ((mode==1) || (collision(pos, shield, 0, NULL)==0)) break;
	  POLY_TRACE(printf("s"); fflush(NULL));
	}
	if (cnt1 >= max_try) { free(poly); return (-1); }
      }
      poly[0] = pos[0]; poly[1] = pos[1]; poly[2] = pos[2];
      max_cnt=imax(cnt1, max_cnt);
      POLY_TRACE(printf("S"); fflush(NULL));
      //POLY_TRACE(/* printf("placed Monomer 0 at (%f,%f,%f)\n",pos[0],pos[1],pos[2]) */);

      poz[0]=pos[0]; poz[1]=pos[1]; poz[2]=pos[2];

      /* place 2nd monomer */
      n=1;
      if (posed2 != NULL && posed != NULL && angle2 > -1.0) {
	/* if position of 2nd monomer is given */
	pos[0]=posed2[0];
	pos[1]=posed2[1];
	pos[2]=posed2[2];
	/* calculate preceding monomer so that bond_length is correct */
	absc=sqrt(SQR(pos[0]-poz[0])+SQR(pos[1]-poz[1])+SQR(pos[2]-poz[2]));
	poz[0]=pos[0]+(poz[0]-pos[0])*bond_length/absc;
	poz[1]=pos[1]+(poz[1]-pos[1])*bond_length/absc;
	poz[2]=pos[2]+(poz[2]-pos[2])*bond_length/absc;
	//POLY_TRACE(/* printf("virtually shifted position of first monomer to (%f,%f,%f)\n",poz[0],poz[1],poz[2]) */);
      } else {
	/* randomly place 2nd monomer */
	for (cnt1=0; cnt1<max_try; cnt1++) {
	  zz     = (2.0*d_random()-1.0)*bond_length;
          rr     = sqrt(SQR(bond_length)-SQR(zz));
	  phi    = 2.0*PI*d_random();
	  pos[0] = poz[0]+rr*cos(phi);
	  pos[1] = poz[1]+rr*sin(phi);
	  pos[2] = poz[2]+zz;
#ifdef CONSTRAINTS
	  if(constr==0 || constraint_collision(pos,poly+3*(n-1))==0){
#endif

	    if (mode==1 || collision(pos, shield, n, poly)==0) break;
	    if (mode==0) { cnt1 = -1; break; }
#ifdef CONSTRAINTS
	  }
#endif
	  POLY_TRACE(printf("m"); fflush(NULL));
	}
	if (cnt1 >= max_try) {
	  fprintf(stderr,"\nWarning! Attempt #%d to build polymer %d failed while placing monomer 2!\n",cnt2+1,p);
	  fprintf(stderr,"         Retrying by re-setting the start-monomer of current chain...\n");
	}
	if (cnt1 == -1 || cnt1 >= max_try) {
	  continue; /* continue the main loop */
	}
      }
      if(posed2!=NULL && p>0) {
	free(posed2);
	posed2=NULL;
      }
      poly[3*n] = pos[0]; poly[3*n+1] = pos[1]; poly[3*n+2] = pos[2];
      max_cnt=imax(cnt1, max_cnt);
      POLY_TRACE(printf("M"); fflush(NULL));
      //POLY_TRACE(/* printf("placed Monomer 1 at (%f,%f,%f)\n",pos[0],pos[1],pos[2]) */);
      
      /* place remaining monomers */
      for (n=2; n<MPC; n++) { 
	if (angle2 > -1.0) {
	  if(n==2) { /* if the 2nd angle is set, construct preceding monomer 
		       with resulting plane perpendicular on the xy-plane */
	    poy[0]=2*poz[0]-pos[0];
	    poy[1]=2*poz[1]-pos[1];
	    if(pos[2]==poz[2])
	      poy[2]=poz[2]+1;
	    else
	      poy[2]=poz[2];
	  } else { 
	    /* save 3rd last monomer */
	    pox[0]=poy[0]; pox[1]=poy[1]; pox[2]=poy[2]; 
	  }
	}
	if (angle > -1.0) { 
	  /* save one but last monomer */
	  poy[0]=poz[0]; poy[1]=poz[1]; poy[2]=poz[2]; 
	}
	/* save last monomer */
	poz[0]=pos[0]; poz[1]=pos[1]; poz[2]=pos[2];

	if(angle > -1.0){
	  a[0]=poy[0]-poz[0];
	  a[1]=poy[1]-poz[1];
	  a[2]=poy[2]-poz[2];

	  b[0]=pox[0]-poy[0];
	  b[1]=pox[1]-poy[1];
	  b[2]=pox[2]-poy[2];

	  vector_product(a,b,c);	  
	}

	for (cnt1=0; cnt1<max_try; cnt1++) {
	  if(angle > -1.0) {
	    if (sqrlen(c) < ROUND_ERROR_PREC) {
	      fprintf(stderr, "WARNING: rotation axis is 0,0,0, check the angles given to the polymer command\n");
	      c[0] = 1; c[1] = 0; c[2] = 0;
	    }
	    if(angle2 > -1.0 && n>2) {
	      vec_rotate(a,angle2,c,d);
	    } else {
	      phi = 2.0*PI*d_random();
	      vec_rotate(a,phi,c,d);
	    }

	    vec_rotate(d,angle,a,b);

	    pos[0] = poz[0] + b[0];
	    pos[1] = poz[1] + b[1];
	    pos[2] = poz[2] + b[2];

	  } else {
            zz     = (2.0*d_random()-1.0)*bond_length;
            rr     = sqrt(SQR(bond_length)-SQR(zz));
            phi    = 2.0*PI*d_random();
            pos[0] = poz[0]+rr*cos(phi);
            pos[1] = poz[1]+rr*sin(phi);
            pos[2] = poz[2]+zz;
	  }
	  
	  //POLY_TRACE(/* printf("a=(%f,%f,%f) absa=%f M=(%f,%f,%f) c=(%f,%f,%f) absMc=%f a*c=%f)\n",a[0],a[1],a[2],sqrt(SQR(a[0])+SQR(a[1])+SQR(a[2])),M[0],M[1],M[2],c[0],c[1],c[2],sqrt(SQR(M[0]+c[0])+SQR(M[1]+c[1])+SQR(M[2]+c[2])),a[0]*c[0]+a[1]*c[1]+a[2]*c[2]) */);
	  //POLY_TRACE(/* printf("placed Monomer %d at (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);

#ifdef CONSTRAINTS
	  if(constr==0 || constraint_collision(pos,poly+3*(n-1))==0){
#endif
	    if (mode==1 || collision(pos, shield, n, poly)==0) break;
	    if (mode==0) { cnt1 = -2; break; }
#ifdef CONSTRAINTS
	  }
#endif
	  POLY_TRACE(printf("m"); fflush(NULL));
	}
	if (cnt1 >= max_try) {
	  fprintf(stderr,"\nWarning! Attempt #%d to build polymer %d failed after %d unsuccessful trials to place monomer %d!\n",cnt2+1,p,cnt1,n);
	  fprintf(stderr,"         Retrying by re-setting the start-monomer of current chain...\n");
	}
	if (cnt1 == -2 || cnt1 >= max_try) {
	  n=0; break;
	}
	poly[3*n] = pos[0]; poly[3*n+1] = pos[1]; poly[3*n+2] = pos[2];
	max_cnt=imax(cnt1, max_cnt);
	POLY_TRACE(printf("M"); fflush(NULL));
      }
      if (n>0) break;
    } /* cnt2 */
    POLY_TRACE(printf(" %d/%d->%d \n",cnt1,cnt2,max_cnt));
    if (cnt2 >= max_try) { free(poly); return(-2); } else max_cnt = imax(max_cnt,imax(cnt1,cnt2));

    /* actually creating current polymer in ESPResSo */
    for (n=0; n<MPC; n++) {
      
      pos[0] = poly[3*n]; pos[1] = poly[3*n+1]; pos[2] = poly[3*n+2];
      if (place_particle(part_id, pos)==ES_PART_ERROR ||
	  (set_particle_q(part_id, ((n % cM_dist==0) ? val_cM : 0.0) )==ES_ERROR) ||
	  (set_particle_type(part_id, ((n % cM_dist==0) ? type_cM : type_nM) )==ES_ERROR))
	{ free(poly); return (-3); }
      
      if(n>=bond_size){
	bond[1] = part_id - bond_size;
	for(i=2;i<=bond_size;i++){
	  bond[i] = part_id - bond_size + i;
	}
	if(change_particle_bond(part_id-bond_size+1, bond, 0)==ES_ERROR)
	  { free(poly); return (-3); }
      }
      part_id++;
      //POLY_TRACE(/* printf("placed Monomer %d at (%f,%f,%f)\n",n,pos[0],pos[1],pos[2]) */);
    }
  }
  free(poly);
  return(imax(max_cnt,cnt2));
}

int counterionsC(int N_CI, int part_id, int mode, double shield, int max_try, double val_CI, int type_CI) {
  int n, cnt1,max_cnt;
  double pos[3];

  cnt1 = max_cnt = 0;
  for (n=0; n<N_CI; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      pos[0]=box_l[0]*d_random();
      pos[1]=box_l[1]*d_random();
      pos[2]=box_l[2]*d_random();
      if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
      POLY_TRACE(printf("c"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==ES_PART_ERROR) return (-3);
    if (set_particle_q(part_id, val_CI)==ES_ERROR) return (-3);
    if (set_particle_type(part_id, type_CI)==ES_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("C"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-1);
  return(imax(max_cnt,cnt1));
}

int saltC(int N_pS, int N_nS, int part_id, int mode, double shield, int max_try, double val_pS, double val_nS, int type_pS, int type_nS, double rad) {
  int n, cnt1,max_cnt;
  double pos[3], dis2;

  cnt1 = max_cnt = 0;

  /* Place positive salt ions */
  for (n=0; n<N_pS; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      if (rad > 0.) {
        pos[0]=rad*(2.*d_random()-1.);
        pos[1]=rad*(2.*d_random()-1.);
        pos[2]=rad*(2.*d_random()-1.);
        dis2 = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
        pos[0] += box_l[0]*0.5;
        pos[1] += box_l[1]*0.5;
        pos[2] += box_l[2]*0.5;
        if (((mode!=0) || (collision(pos, shield, 0, NULL)==0)) && (dis2 < (rad * rad))) break;
      } else {
        pos[0]=box_l[0]*d_random();
        pos[1]=box_l[1]*d_random();
        pos[2]=box_l[2]*d_random();
        if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
      }
      POLY_TRACE(printf("p"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==ES_PART_ERROR) return (-3);
    if (set_particle_q(part_id, val_pS)==ES_ERROR) return (-3);
    if (set_particle_type(part_id, type_pS)==ES_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("P"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-1);

  /* Place negative salt ions */
  for (n=0; n<N_nS; n++) {
    for (cnt1=0; cnt1<max_try; cnt1++) {
      if (rad > 0.) {
        pos[0]=rad*(2.*d_random()-1.);
        pos[1]=rad*(2.*d_random()-1.);
        pos[2]=rad*(2.*d_random()-1.);
        dis2 = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
        pos[0] += box_l[0]*0.5;
        pos[1] += box_l[1]*0.5;
        pos[2] += box_l[2]*0.5;
        if (((mode!=0) || (collision(pos, shield, 0, NULL)==0)) && (dis2 < (rad * rad))) break;
      } else {
        pos[0]=box_l[0]*d_random();
        pos[1]=box_l[1]*d_random();
        pos[2]=box_l[2]*d_random();
        if ((mode!=0) || (collision(pos, shield, 0, NULL)==0)) break;
      }
      POLY_TRACE(printf("n"); fflush(NULL));
    }
    if (cnt1 >= max_try) return (-1);
    if (place_particle(part_id, pos)==ES_PART_ERROR) return (-3);
    if (set_particle_q(part_id, val_nS)==ES_ERROR) return (-3);
    if (set_particle_type(part_id, type_nS)==ES_ERROR) return (-3);
    part_id++; max_cnt=imax(cnt1, max_cnt);
    POLY_TRACE(printf("N"); fflush(NULL));
  }
  POLY_TRACE(printf(" %d->%d \n",cnt1,max_cnt));
  if (cnt1 >= max_try) return(-2);
  return(imax(max_cnt,cnt1));
}

double velocitiesC(double v_max, int part_id, int N_T) {
  double v[3], v_av[3];
  int i;

  v_av[0] = v_av[1] = v_av[2] = 0.0;
  for (i=part_id; i < part_id+N_T; i++) {
    do {
      v[0] = v_max * 2.*(d_random()-.5) * time_step;
      v[1] = v_max * 2.*(d_random()-.5) * time_step;
      v[2] = v_max * 2.*(d_random()-.5) * time_step;
      // note that time_step == -1, as long as it is not yet set
    } while ( sqrt(SQR(v[0])+SQR(v[1])+SQR(v[2])) > v_max * fabs(time_step));
    v_av[0]+=v[0]; v_av[1]+=v[1]; v_av[2]+=v[2];
    if (set_particle_v(i, v)==ES_ERROR) {
      fprintf(stderr, "INTERNAL ERROR: failed upon setting one of the velocities in Espresso (current average: %f)!\n",
	      sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2]))); 
      fprintf(stderr, "Aborting...\n"); errexit();
    }
  }
  // note that time_step == -1, as long as it is not yet set
  return ( sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2])) / fabs(time_step) );
}

double maxwell_velocitiesC(int part_id, int N_T) {
  double v[3], v_av[3],uniran[2];
  int i;
  int flag=1;
  uniran[0]=d_random();
  uniran[1]=d_random();
  v_av[0] = v_av[1] = v_av[2] = 0.0;
  for (i=part_id; i < part_id+N_T; i++) {
    if(flag == 1 ) {
      v[0] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      v[1] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      uniran[0]=d_random();
      uniran[1]=d_random();
      v[2] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      flag = 0;
    } else {
      v[0] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      uniran[0]=d_random();
      uniran[1]=d_random();
      v[1] = pow((-2. * log(uniran[0])),0.5) * cos (2. * PI * uniran[1]) * time_step;
      v[2] = pow((-2. * log(uniran[1])),0.5) * sin (2. * PI * uniran[0]) * time_step;
      flag = 1;      
    }
    //printf("%f \n %f \n %f \n",v[0],v[1],v[2]);
    v_av[0]+=v[0]; v_av[1]+=v[1]; v_av[2]+=v[2];
    if (set_particle_v(i, v)==ES_ERROR) {
      fprintf(stderr, "INTERNAL ERROR: failed upon setting one of the velocities in Espresso (current average: %f)!\n",sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2]))); 
      fprintf(stderr, "Aborting...\n"); errexit();
    }
  }
  // note that time_step == -1, as long as it is not yet set
  return ( sqrt(SQR(v_av[0])+SQR(v_av[1])+SQR(v_av[2])) / fabs(time_step) );
}

int collectBonds(int mode, int part_id, int N_P, int MPC, int type_bond, int **bond_out, int ***bonds_out) {
  int i,j,k,ii,size, *bond=NULL, **bonds=NULL;

  /* Get particle and bonding informations. */
  IntList *bl;
  Particle *prt, *sorted;
  bl  = (IntList*)malloc(1*sizeof(IntList));
  prt = (Particle*)malloc(n_part*sizeof(Particle));
  mpi_get_particles(prt, bl); 

  /* Sort the received informations. */
  sorted = (Particle*)malloc(n_part*sizeof(Particle));
  for(i = 0; i < n_part; i++)
    memcpy(&sorted[prt[i].p.identity], &prt[i], sizeof(Particle));
  free(prt);
  prt = sorted;
  
  if (mode == 1) {
    /* Find all the bonds leading to and from the ending monomers of the chains. */
    bond  = (int*)malloc(2*N_P*sizeof(int));      
    bonds   = (int**)malloc(2*N_P*sizeof(int *));
    for (i=0; i < 2*N_P; i++) { 
      bond[i]=0;  
      bonds[i]= (int*)malloc(1*sizeof(int)); 
    }
    for (k=part_id; k < N_P*MPC + part_id; k++) {
      i=0;
      while(i < prt[k].bl.n) {
	size = bonded_ia_params[prt[k].bl.e[i]].num;
	if (prt[k].bl.e[i++] == type_bond) {
	  for(j=0; j<size; j++) {
	    if ((prt[k].p.identity % MPC == 0) || ( (prt[k].p.identity+1) % MPC == 0)) {
	      ii = prt[k].p.identity%MPC ? 2*(prt[k].p.identity+1)/MPC-1 : 2*prt[k].p.identity/MPC;
	      bonds[i] = (int*)realloc(bonds[i], (bond[i]+1)*sizeof(int));
	      bonds[ii][bond[ii]++] = prt[k].bl.e[i];
	    }
	    else if ((prt[k].bl.e[i] % MPC == 0) || ( (prt[k].bl.e[i]+1) % MPC == 0)) {
	      ii = prt[k].bl.e[i]%MPC ? 2*(prt[k].bl.e[i]+1)/MPC-1 : 2*prt[k].bl.e[i]/MPC;
	      bonds[i] = (int*)realloc(bonds[i], (bond[i]+1)*sizeof(int));
	      bonds[ii][bond[ii]++] = prt[k].p.identity;
	    }
	    i++;
	  }
	}
	else i += size;
      }
    }
    POLY_TRACE(for (i=0; i < 2*N_P; i++) {
      printf("(%d) %d:\t",i,i%2 ? (i+1)*MPC/2-1 : i*MPC/2); if(bond[i]>0) for(j=0;j<bond[i];j++) printf("%d ",bonds[i][j]); printf("\t=%d\n",bond[i]);
    });
  }
  else if (mode == 2) {
    /* Find all the bonds leading to and from each monomer. */
    bond  = (int*)malloc(N_P*MPC*sizeof(int));                
    bonds   = (int**)malloc(N_P*MPC*sizeof(int *));
    for (i=0; i < N_P*MPC + part_id; i++) { 
      bond[i]=0;  
      bonds[i]= (int*)malloc(1*sizeof(int)); 
    }
    for (k=part_id; k < N_P*MPC + part_id; k++) {
      i=0;
      while(i < prt[k].bl.n) {
	size = bonded_ia_params[prt[k].bl.e[i]].num;
	if (prt[k].bl.e[i++] == type_bond) {
	  for(j=0; j<size; j++) {
	    ii = prt[k].bl.e[i];
	    bonds[k] = (int*) realloc(bonds[k], (bond[k]+1)*sizeof(int));
	    bonds[k][bond[k]++] = ii;
	    bonds[ii] = (int*) realloc(bonds[ii], (bond[ii]+1)*sizeof(int));
	    bonds[ii][bond[ii]++] = k;
	    i++;
	  }
	}
	else i += size;
      }
    }
    POLY_TRACE(for (i=0; i < N_P*MPC + part_id; i++) { 
      printf("%d:\t",i); if(bond[i]>0) for(j=0;j<bond[i];j++) printf("%d ",bonds[i][j]); printf("\t=%d\n",bond[i]); 
    });
  }
  else {
    fprintf(stderr, "Unknown mode %d requested!\nAborting...\n",mode); fflush(NULL); return(-2);
  }
  free(prt); realloc_intlist(bl, 0);
  *bond_out  = bond;
  *bonds_out = bonds;
  return(0);
}

int crosslinkC(int N_P, int MPC, int part_id, double r_catch, int link_dist, int chain_dist, int type_bond, int max_try) {
  int i,j,k,ii,size, bondN[2], *bond, **bonds, *link, **links, *cross, crossL;

  /* Find all the bonds leading to and from each monomer. */
  if (collectBonds(2, part_id, N_P, MPC, type_bond, &bond, &bonds)) return(-2);
  POLY_TRACE(for (i=0; i < N_P*MPC + part_id; i++) { 
      printf("%d:\t",i); 
      if(bond[i]>0) 
        for(j=0;j<bond[i];j++) 
          printf("%d ",bonds[i][j]); 
      printf("\t=%d\n",bond[i]); 
    });
  
  /* Find all possible binding partners in the neighbourhood of the unconnected ending monomers. */
  link  = (int*)malloc(2*N_P*sizeof(int));       
  links = (int**)malloc(2*N_P*sizeof(int *));
  for (i=0; i < N_P; i++) {
    for (k=0; k<2; k++) {
      if (bond[i*MPC+k*(MPC-1)] == 1) {
	links[2*i+k] = (int*)malloc(n_part*sizeof(int));
	link[2*i+k] = mindist3(i*MPC+k*(MPC-1)+part_id, r_catch, links[2*i+k]);
	links[2*i+k] = (int*)realloc(links[2*i+k],link[2*i+k]*sizeof(int));
      }
      else if (bond[i*MPC+k*(MPC-1)] == 2) link[2*i+k] = -1;  /* Note that links[2*i+k] will not be malloc()ed now (taken care of at end)!!! */
      else { fprintf(stderr,"Runaway end-monomer %d detected (has %d bonds)!\nAborting...\n", i*N_P+k*(MPC-1)+part_id, bond[i*MPC+k*(MPC-1)]); 
             fflush(NULL); return(-2); }
      POLY_TRACE(printf("%d: ",i*MPC+k*(MPC-1)+part_id); 
		 for (j=0; j<link[2*i+k]; j++) printf("%d ",links[2*i+k][j]); printf("\t=%d\n",link[2*i+k]); fflush(NULL) );
    }
  }

  /* Throw out all the monomers which are ends, which are too close to the ending monomers on the same chain, or which are no monomers at all. */
  for (i=0; i < N_P; i++) {
    for (k=0; k<2; k++) {
      size = 0;  ii = i*MPC + k*(MPC-1) + part_id;
      if (link[2*i+k] >= 0) {
	for (j=0; j < link[2*i+k]; j++) {                     /* only monomers && ((same chain, but sufficiently far away) || (different chain)) */
	  if ( (links[2*i+k][j] < N_P*MPC+part_id) && ( ((abs(links[2*i+k][j] - ii) > chain_dist) || (abs(links[2*i+k][j]-i*MPC) > (1.*MPC))) ) )
	    if ((links[2*i+k][j] % MPC != 0) && ((links[2*i+k][j]+1) % MPC != 0)) links[2*i+k][size++] = links[2*i+k][j];    /* no ends accepted */
	}
	link[2*i+k]  = size; 
	links[2*i+k] = (int*)realloc(links[2*i+k],link[2*i+k]*sizeof(int));
      }
      POLY_TRACE(printf("%d: ",ii); for (j=0; j<link[2*i+k]; j++) printf("%d ",links[2*i+k][j]); printf("\t=%d\n",link[2*i+k]); fflush(NULL) );
    }
  }

  /* Randomly choose a partner (if not available -> '-1') for each polymer chain's end if it's not already been crosslinked (-> '-2'). */
  cross = (int*)malloc(2*N_P*sizeof(int)); crossL = 0;
  for (i=0; i < 2*N_P; i++) 
    if (link[i] > 0) { cross[i] = links[i][(int)dround(d_random()*(link[i]-1))]; crossL++; }  else { cross[i] = -1+link[i]; crossL -= link[i]; }
  POLY_TRACE(for (i=0; i < 2*N_P; i++) printf("%d -> %d \t", i%2 ? (i+1)*MPC/2-1 : i*MPC/2, cross[i]); printf("=> %d\n",crossL); fflush(NULL) );

  /* Remove partners (-> '-3') if they are less than link_dist apart and retry. */
  k = 0; ii = 1;
  while ((k < max_try) && (ii > 0)) {
    POLY_TRACE(printf("Check #%d: ",k));
    for (i=0; i < 2*N_P; i++) {
      if (cross[i] >= 0) {
	for (j=0; j < 2*N_P; j++) {       /* In the neighbourhood of each partner shall be no future crosslinks (preventing stiffness). */
	  if ((j != i) && (cross[j] >=0) && (abs(cross[j]-cross[i]) < link_dist)) {
	    cross[i] = -3; cross[j] = -3; crossL -= 2; POLY_TRACE(printf("%d->%d! ",i,j)); break;
	  }
	}
	if (cross[i] == -3) continue;     /* Partners shall not be too close to the chain's ends (because these will be crosslinked at some point). */
	if ((cross[i] % MPC < link_dist) || (cross[i] % MPC >= MPC-link_dist)) {
	  cross[i] = -3; crossL--; POLY_TRACE(printf("%d->end! ",i)); }      
	else {                            /* In the neighbourhood of each partner there shall be no other crosslinks (preventing stiffness). */
	  for (j = cross[i]-link_dist+1; j < cross[i]+link_dist-1; j++) {
	    if ((j % MPC == 0) || ((j+1) % MPC == 0)) size = 1; else size = 2;
	    if ((bond[j] > size) && (j - floor(i/2.)*MPC < MPC)) {
	      cross[i] = -3; crossL--; POLY_TRACE(printf("%d->link! ",i)); break; 
	    }
	  }
	}
      }
    }   POLY_TRACE(printf("complete => %d CL left; ",crossL));
    if (k == max_try-1) break; else ii = 0;  /* Get out if max_try is about to be reached, preventing dangling unchecked bond suggestions. */
    if (crossL < 2*N_P) {
      for (i=0; i < 2*N_P; i++) {         /* If crosslinks violated the rules & had to be removed, create new ones now. */
	if (cross[i] == -3) {
	  ii++;
	  if (link[i] > 0) { cross[i] = links[i][(int)dround(d_random()*(link[i]-1))]; crossL++; } else { return(-2); }
	}
      }
    }   POLY_TRACE(printf("+ %d new = %d CL.\n",ii,crossL));
    if (ii > 0) k++;
  }
  POLY_TRACE(for (i=0; i < 2*N_P; i++) printf("%d -> %d \t", i%2 ? (i+1)*MPC/2-1 : i*MPC/2, cross[i]); printf("=> %d\n",crossL); fflush(NULL) );

  /* Submit all lawful partners as new bonds to Espresso (observing that bonds are stored with the higher-ID particle only). */
  if (k >= max_try) return(-1);

  {
    size = 0;
    for (i=0; i < N_P; i++) {
      if (cross[2*i] >= 0 ) {
	bondN[0] = type_bond; bondN[1] = i*MPC + part_id; size++;
	if (change_particle_bond(cross[2*i], bondN, 0)==ES_ERROR) return (-3);
      }
      if (cross[2*i+1] >= 0) {
	bondN[0] = type_bond; bondN[1] = cross[2*i+1]; size++;
	if (change_particle_bond(i*MPC+(MPC-1) + part_id, bondN, 0)==ES_ERROR) return (-3);
      }
      free(bonds[2*i]);    if (link[2*i]   >= 0) free(links[2*i]);    /* else crash(); because links[2*i]   has never been malloc()ed then */
      free(bonds[2*i+1]);  if (link[2*i+1] >= 0) free(links[2*i+1]);  /* else crash(); because links[2*i+1] has never been malloc()ed then */
    }
    free(bond); free(bonds); free(link); free(links); free(cross);
    POLY_TRACE(printf("Created %d new bonds; now %d ends are crosslinked!\n", size, crossL));
    return(crossL); 
  }
}

int diamondC(double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet) {
  int i,j,k, part_id, bond[2], type_bond=0,type_node=0,type_cM=1,type_nM=1, type_CI=2;
  double pos[3], off = bond_length/sqrt(3);
  double dnodes[8][3]  = {{0,0,0}, {1,1,1}, {2,2,0}, {0,2,2}, {2,0,2}, {3,3,1}, {1,3,3}, {3,1,3}};
  int    dchain[16][5] = {{0,1, +1,+1,+1}, {1,2, +1,+1,-1}, {1,3, -1,+1,+1}, {1,4, +1,-1,+1},
			  {2,5, +1,+1,+1}, {3,6, +1,+1,+1}, {4,7, +1,+1,+1}, {5,0, +1,+1,-1},
			  {5,3, +1,-1,+1}, {5,4, -1,+1,+1}, {6,0, -1,+1,+1}, {6,2, +1,-1,+1},
			  {6,4, +1,+1,-1}, {7,0, +1,-1,+1}, {7,2, -1,+1,+1}, {7,3, +1,+1,-1}};

  part_id = 0;
  /* place 8 tetra-functional nodes */
  for(i=0; i<8; i++) {
    for(j=0; j<3; j++) { 
      dnodes[i][j] *= a/4.; pos[j] = dnodes[i][j]; 
    }
    if (place_particle(part_id, pos)==ES_PART_ERROR) return (-3);
    if (set_particle_q(part_id, val_nodes)==ES_ERROR) return (-3);
    if (set_particle_type(part_id, type_node)==ES_ERROR) return (-3);
    part_id++;
  }

  /* place intermediate monomers on chains connecting the nodes */
  for(i=0; i<2*8; i++) {
    for(k=1; k<=MPC; k++) {
      for(j=0; j<3; j++) pos[j] = dnodes[dchain[i][0]][j] + k*dchain[i][2+j]*off;
      if (place_particle(part_id, pos)==ES_PART_ERROR) return (-3);
      if (set_particle_q(part_id, (k % cM_dist==0) ? val_cM : 0.0)==ES_ERROR) return (-3);
      if (set_particle_type(part_id, (k % cM_dist==0) ? type_cM : type_nM)==ES_ERROR) return (-3);
      bond[0] = type_bond; 
      if(k==1) { 
	if(nonet!=1) { bond[1] = dchain[i][0]; if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3); } }
      else { 
	bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3); }
      if((k==MPC)&&(nonet!=1)) { 
	bond[1] = dchain[i][1];
	if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3);
      }
      part_id++;
    }
  }

  /* place counterions (if any) */
  if(N_CI > 0) counterionsC(N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);
  
  return(0);
}

int icosaederC(double ico_a, int MPC, int N_CI, double val_cM, double val_CI, int cM_dist) {
  int i,j,k,l, part_id, bond[2], type_bond=0,type_cM=0,type_nM=1, type_CI=2;
  double pos[3],pos_shift[3], vec[3],e_vec[3],vec_l, bond_length=(2*ico_a/3.)/(1.*MPC);
  double ico_g=ico_a*(1+sqrt(5))/2.0, shift=0.0;
  double ico_coord[12][3] = {{0,+ico_a,+ico_g}, {0,+ico_a,-ico_g}, {0,-ico_a,+ico_g}, {0,-ico_a,-ico_g},
			     {+ico_a,+ico_g,0}, {+ico_a,-ico_g,0}, {-ico_a,+ico_g,0}, {-ico_a,-ico_g,0},
			     {+ico_g,0,+ico_a}, {-ico_g,0,+ico_a}, {+ico_g,0,-ico_a}, {-ico_g,0,-ico_a}};
  int    ico_NN[12][5]    = {{2,8,4,6, 9}, {3,10,4,6,11}, {0,8,5,7, 9}, {1,10,5,7,11}, {0,6,1,10,8}, {2,7,3,10,8},
			     {0,4,1,11,9}, {2,5,3,11, 9}, {0,2,5,10,4}, {0,2,7,11, 6}, {1,3,5,8, 4}, {1,3,7,9, 6}};
  int    ico_ind[12][10];

  /* make sure that the edges in ico_NN are sorted such that NearestNeighbours are next to each other */
  /* int    ico_NN[12][5]    = {{2,4,6,8, 9}, {3,4,6,10,11}, {0,5,7,8, 9}, {1,5,7,10,11}, {0,1,6,8,10}, {2,3,7,8,10}, 
			     {0,1,4,9,11}, {2,3,5, 9,11}, {0,2,4,5,10}, {0,2,6, 7,11}, {1,3,4,5, 8}, {1,3,6,7, 9}};
  for(i=0; i<12; i++) {
    printf("%d: { ",i); 
    for(j=0; j<5; j++) printf("%d ",ico_NN[i][j]);
    printf("} -> ");
    for(j=0; j<5; j++) 
      for(l=0; l<5; l++) 
	if(j!=l) 
	  for(k=0; k<5; k++) 
	    if(ico_NN[ico_NN[i][j]][k]==ico_NN[i][l]) printf("%d = %d (@%d)  ",ico_NN[i][j],ico_NN[i][l],k);
    printf("\n");
  } */

  /* shift coordinates to not be centered around zero but rather be positive */
  if(ico_a > ico_g) shift=ico_a; else shift=ico_g;

  /* create fulleren & soccer-ball */
  part_id = 0;
  for(i=0; i<12; i++) {
    for(j=0; j<5; j++) {
      /* place chains along the 5 edges around each of the 12 icosaeder's vertices */
      if(j < 4) for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][j+1]][l] - ico_coord[ico_NN[i][j]][l])/3.;
      else      for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][0]][l]   - ico_coord[ico_NN[i][4]][l])/3.;
      vec_l = sqrt(SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]));
      for(l=0; l<3; l++) e_vec[l] = vec[l]/vec_l;

      ico_ind[i][j] = part_id; bond_length = vec_l/(1.*MPC);
      for(l=0; l<3; l++) pos[l] = ico_coord[i][l] + (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3.;
      for(k=0; k<MPC; k++) {
	for(l=0; l<3; l++) pos_shift[l] = pos[l] + shift;
	if (place_particle(part_id, pos_shift)==ES_PART_ERROR) return (-3);
	if (set_particle_q(part_id, val_cM)==ES_ERROR) return (-3);
	if (set_particle_type(part_id, type_cM)==ES_ERROR) return (-3);
	bond[0] = type_bond;
	if (k > 0) {
	  bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3); 
	}
	part_id++;
	for(l=0; l<3; l++) pos[l] += bond_length*e_vec[l];
      }

      /* place chains along the 5 edges on the middle third of the connection between two NN vertices */
      if(i < ico_NN[i][j]) {
	for(l=0; l<3; l++) vec[l] = (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3.;
	vec_l = sqrt(SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]));
	for(l=0; l<3; l++) e_vec[l] = vec[l]/vec_l;
	
	ico_ind[i][j+5] = part_id; bond_length = vec_l/(1.*MPC);
	for(l=0; l<3; l++) pos[l] = ico_coord[i][l] + (ico_coord[ico_NN[i][j]][l] - ico_coord[i][l])/3. + bond_length*e_vec[l];
	for(k=1; k<MPC; k++) {
	  for(l=0; l<3; l++) pos_shift[l] = pos[l] + shift;
	  if (place_particle(part_id, pos_shift)==ES_ERROR) return (-3);
	  if (set_particle_q(part_id, 0.0)==ES_ERROR) return (-3);
	  if (set_particle_type(part_id, type_nM)==ES_ERROR) return (-3);
	  bond[0] = type_bond;
	  if (k > 1) {
	    bond[1] = part_id-1; if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3); }
	  else {
	    bond[1] = ico_ind[i][j]; if (change_particle_bond(part_id, bond, 0)==ES_ERROR) return (-3); }
	  part_id++;
	  for(l=0; l<3; l++) pos[l] += bond_length*e_vec[l];
	}
      } 
    }

    for(j=0; j<5; j++) {
      /* add bonds between the edges around the vertices */
      bond[0] = type_bond;
      //      if(j>0) bond[1] = ico_ind[i][j-1] + (MPC-1); else bond[1] = ico_ind[i][4] + (MPC-1);
      if(j>0) bond[1] = ico_ind[i][j-1] + (MPC-1); else if(MPC>0) bond[1] = ico_ind[i][4] + (MPC-1); else bond[1] = ico_ind[i][4];
      if (change_particle_bond(ico_ind[i][j], bond, 0)==ES_ERROR) return (-2);

      /* connect loose edges around vertices with chains along the middle third already created earlier */
      if(i > ico_NN[i][j]) {
	bond[0] = type_bond;
	for(l=0; l<5; l++) if(ico_NN[ico_NN[i][j]][l] == i) break;
	if(l==5) {
	  fprintf(stderr, "INTERNAL ERROR: Couldn't find my neighbouring edge upon creating the icosaeder!\n");
	  errexit();
	}
	bond[1] = ico_ind[ico_NN[i][j]][l+5] + (MPC-2);
	if (change_particle_bond(ico_ind[i][j], bond, 0)==ES_ERROR) return (-1);
      }
    }
  }

  /* place counterions (if any) */
  if(N_CI > 0) counterionsC(N_CI, part_id, 1, 0.0, 30000, val_CI, type_CI);
  
  return(0);
}

