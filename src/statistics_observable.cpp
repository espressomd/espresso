/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  
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
#include "statistics_observable.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"
#include "lb.hpp"
#include "pressure.hpp"
#include "rotation.hpp"

observable** observables = 0;
int n_observables = 0; 

int observable_particle_velocities(void* idlist, double* A, unsigned int n_A) {
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

int observable_particle_angular_momentum(void* idlist, double* A, unsigned int n_A) {
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for ( int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;

    #ifdef ROTATION

    double RMat[9];
    double omega[3];
    define_rotation_matrix(&partCfg[ids->e[i]], RMat);
    omega[0] = RMat[0 + 3*0]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*0]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*0]*partCfg[ids->e[i]].m.omega[2];
    omega[1] = RMat[0 + 3*1]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*1]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*1]*partCfg[ids->e[i]].m.omega[2];
    omega[2] = RMat[0 + 3*2]*partCfg[ids->e[i]].m.omega[0] + RMat[1 + 3*2]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*2]*partCfg[ids->e[i]].m.omega[2];

    A[3*i + 0] = omega[0];
    A[3*i + 1] = omega[1];
    A[3*i + 2] = omega[2];

    #else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

    #endif
  }
  return 0;
}

#ifdef ELECTROSTATICS
int observable_particle_currents(void* idlist, double* A, unsigned int n_A) {
  double charge;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    A[3*i + 0] = charge * partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = charge * partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

int observable_currents(void* idlist, double* A, unsigned int n_A) {
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_part)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].m.v[0]/time_step;
    j[1] += charge * partCfg[ids->e[i]].m.v[1]/time_step;
    j[2] += charge * partCfg[ids->e[i]].m.v[2]/time_step;
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}

int observable_dipole_moment(void* idlist, double* A, unsigned int n_A) {
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] > n_part)
      return 1;
    charge = partCfg[ids->e[i]].p.q;
    j[0] += charge * partCfg[ids->e[i]].r.p[0];
    j[1] += charge * partCfg[ids->e[i]].r.p[1];
    j[2] += charge * partCfg[ids->e[i]].r.p[2];
  }
  A[0]=j[0];
  A[1]=j[1];
  A[2]=j[2];
  return 0;
}
#endif

int observable_com_velocity(void* idlist, double* A, unsigned int n_A) {
  double v_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    v_com[0] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].m.v[0]/time_step;
    v_com[1] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].m.v[1]/time_step;
    v_com[2] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].m.v[2]/time_step;
    total_mass += PMASS(partCfg[ids->e[i]]);
  }
  A[0]=v_com[0]/total_mass;
  A[1]=v_com[1]/total_mass;
  A[2]=v_com[2]/total_mass;
  return 0;
}

int observable_blocked_com_velocity(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  n_blocks=n_A/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    total_mass = 0;
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_part)
        return 1;
      A[3*block+0] +=  PMASS(partCfg[id])*partCfg[id].m.v[0]/time_step;
      A[3*block+1] +=  PMASS(partCfg[id])*partCfg[id].m.v[1]/time_step;
      A[3*block+2] +=  PMASS(partCfg[id])*partCfg[id].m.v[2]/time_step;
      total_mass += PMASS(partCfg[ids->e[i]]);
    }
    A[3*block+0] /=  total_mass;
    A[3*block+1] /=  total_mass;
    A[3*block+2] /=  total_mass;
  }
  return 0;
}

int observable_blocked_com_position(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  n_blocks=n_A/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    total_mass = 0;
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_part)
        return 1;
      A[3*block+0] +=  PMASS(partCfg[id])*partCfg[id].r.p[0];
      A[3*block+1] +=  PMASS(partCfg[id])*partCfg[id].r.p[1];
      A[3*block+2] +=  PMASS(partCfg[id])*partCfg[id].r.p[2];
      total_mass += PMASS(partCfg[ids->e[i]]);
    }
    A[3*block+0] /=  total_mass;
    A[3*block+1] /=  total_mass;
    A[3*block+2] /=  total_mass;
  }
  return 0;
}

int observable_com_position(void* idlist, double* A, unsigned int n_A) {
  double p_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    p_com[0] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].r.p[0];
    p_com[1] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].r.p[1];
    p_com[2] += PMASS(partCfg[ids->e[i]])*partCfg[ids->e[i]].r.p[2];
    total_mass += PMASS(partCfg[ids->e[i]]);
  }
  A[0]=p_com[0]/total_mass;
  A[1]=p_com[1]/total_mass;
  A[2]=p_com[2]/total_mass;
  return 0;
}


int observable_com_force(void* idlist, double* A, unsigned int n_A) {
  double f_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    f_com[0] += partCfg[ids->e[i]].f.f[0]/time_step/time_step*2;
    f_com[1] += partCfg[ids->e[i]].f.f[1]/time_step/time_step*2;
    f_com[2] += partCfg[ids->e[i]].f.f[2]/time_step/time_step*2;
  }
  A[0]=f_com[0];
  A[1]=f_com[1];
  A[2]=f_com[2];
  return 0;
}


int observable_blocked_com_force(void* idlist, double* A, unsigned int n_A) {
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  n_blocks=n_A/3; 
  blocksize=ids->n/n_blocks;
  for ( block = 0; block < n_blocks; block++ ) {
    for ( i = 0; i < blocksize; i++ ) {
      id = ids->e[block*blocksize+i];
      if (ids->e[i] >= n_part)
        return 1;
      A[3*block+0] +=  partCfg[id].f.f[0]/time_step/time_step*2;
      A[3*block+1] +=  partCfg[id].f.f[1]/time_step/time_step*2;
      A[3*block+2] +=  partCfg[id].f.f[2]/time_step/time_step*2;
    }
  }
  return 0;
}


int observable_density_profile(void* pdata_, double* A, unsigned int n_A) {
  int binx, biny, binz;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
  double bin_volume=(pdata->maxx-pdata->minx)*(pdata->maxy-pdata->miny)*(pdata->maxz-pdata->minz)/pdata->xbins/pdata->ybins/pdata->zbins;
    
  for (unsigned i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    binx= (int) floor( pdata->xbins*  (ppos[0]-pdata->minx)/(pdata->maxx-pdata->minx));
    biny= (int) floor( pdata->ybins*  (ppos[1]-pdata->miny)/(pdata->maxy-pdata->miny));
    binz= (int) floor( pdata->zbins*  (ppos[2]-pdata->minz)/(pdata->maxz-pdata->minz));
    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      A[binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz] += 1./bin_volume;
    } 
  }
  return 0;
}

#ifdef LB
int observable_lb_velocity_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double xoffset, yoffset, zoffset;
  double x_incr, y_incr, z_incr;
  double p[3], v[3];
  profile_data* pdata;
  pdata=(profile_data*) pdata_;
  int linear_index;

    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  double normalization_factor = 1.;
  if ( pdata->xbins == 1 ) {
    maxi = (int) floor(box_l[0]/lbpar.agrid);
    normalization_factor/=maxi;
    xoffset=0;
    x_incr=lbpar.agrid;
  } else {
    maxi = pdata->xbins;
    xoffset=pdata->minx;
    x_incr=(pdata->maxx-pdata->minx)/(pdata->xbins-1);
  }
  if ( pdata->ybins == 1 ) {
    maxj = (int) floor(box_l[1]/lbpar.agrid);
    normalization_factor/=maxj;
    yoffset=0;
    y_incr=lbpar.agrid;
  } else {
    maxj = pdata->ybins;
    yoffset=pdata->miny;
    y_incr=(pdata->maxy-pdata->miny)/(pdata->ybins-1);
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(box_l[2]/lbpar.agrid);
    normalization_factor/=maxk;
    zoffset=0;
    z_incr=lbpar.agrid;
  } else {
    maxk = pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        p[0]=xoffset + i*x_incr;
        p[1]=yoffset + j*y_incr;
        p[2]=zoffset + k*z_incr;
        if (lb_lbfluid_get_interpolated_velocity(p, v)!=0)
          return 1;
        linear_index = 0;
        if (pdata->xbins > 1)
          linear_index += i*pdata->ybins*pdata->zbins;
        if (pdata->ybins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;

        A[3*linear_index+0]+=v[0];
        A[3*linear_index+1]+=v[1];
        A[3*linear_index+2]+=v[2];
      }
    }
  }
  
  for ( i = 0; i<n_A; i++ ) {
    A[i]*=normalization_factor;
  }

  
  return 0;
}
#endif

#ifdef LB
int observable_lb_radial_velocity_profile(void* pdata_, double* A, unsigned int n_A) {
  unsigned int i, j, k;
  unsigned int maxi, maxj, maxk;
  double roffset, phioffset, zoffset;
  double r, phi, z;
  double r_incr, phi_incr, z_incr;
  double p[3], v[3];
  double v_r, v_phi, v_z;
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  int linear_index;

    
  for ( i = 0; i<n_A; i++ ) {
    A[i]=0;
  }
  double normalization_factor = 1.;
  if ( pdata->rbins == 1 ) {
    return 1;
  } else {
    maxi = pdata->rbins;
    roffset=pdata->minr;
    r_incr=(pdata->maxr-pdata->minr)/(pdata->rbins-1);
  }
  if ( pdata->phibins == 1 ) {
    maxj = (int)floor( 2*3.1415*pdata->maxr/lbpar.agrid ) ; 
    normalization_factor/=maxj;
    phioffset=0;
    phi_incr=2*3.1415/maxj;
  } else {
    maxj = pdata->phibins;
    phioffset=pdata->minphi;
    phi_incr=(pdata->maxphi-pdata->minphi)/(pdata->phibins-1);
  }
  if ( pdata->zbins == 1 ) {
    maxk = (int) floor(box_l[2]/lbpar.agrid);
    normalization_factor/=maxk;
    zoffset=-pdata->center[2];
    z_incr=lbpar.agrid;
  } else {
    maxk = pdata->zbins;
    zoffset=pdata->minz;
    z_incr=(pdata->maxz-pdata->minz)/(pdata->zbins-1);
  }

  for ( i = 0; i < maxi; i++ ) {
    for ( j = 0; j < maxj; j++ ) {
      for ( k = 0; k < maxk; k++ ) {
        r = roffset + i*r_incr;
        phi = phioffset + j*phi_incr;
        z = zoffset + k*z_incr;
        p[0]=r*cos(phi)+pdata->center[0];
        p[1]=r*sin(phi)+pdata->center[1];
        p[2]=z+pdata->center[2];
        if (lb_lbfluid_get_interpolated_velocity(p, v)!=0)
          return 1;
        linear_index = 0;
        if (pdata->rbins > 1)
          linear_index += i*pdata->phibins*pdata->zbins;
        if (pdata->phibins > 1)
          linear_index += j*pdata->zbins;
        if (pdata->zbins > 1)
          linear_index +=k;
        if (r>0) {
          v_r = 1/r*((p[0]-pdata->center[0])*v[0] + (p[1]-pdata->center[1])*v[1]); 
          v_phi = 1/r/r*((p[0]-pdata->center[0])*v[1]-(p[1]-pdata->center[1])*v[0]);
        } else {
          v_r = 0;
          v_phi = 0;
        }
        v_z = v[2];

        A[3*linear_index+0]+=v_r;
        A[3*linear_index+1]+=v_phi;
        A[3*linear_index+2]+=v_z;
      }
    }
  }
  
  for ( i = 0; i<n_A; i++ ) {
    A[i]*=normalization_factor;
  }

  
  return 0;
}
#endif

void transform_to_cylinder_coordinates(double x, double y, double z_, double* r, double* phi, double* z) {
  *z =  z_;
  *r =  sqrt(x*x+y*y);
  *phi = atan2(y,x);
}

int observable_radial_density_profile(void* pdata_, double* A, unsigned int n_A) {
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
    
  for (unsigned i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    transform_to_cylinder_coordinates(ppos[0]-pdata->center[0], ppos[1]-pdata->center[1], ppos[2]-pdata->center[2], &r, &phi, &z);
    //printf("%f %f %f %f %f %f\n", ppos[0], ppos[1], ppos[2], r*cos(phi)+pdata->center[0], r*sin(phi)+pdata->center[1], z+pdata->center[2]);
    binr  =(int)floor((r-pdata->minr)/rbinsize);
    binphi=(int)floor((phi-pdata->minphi)/phibinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);

    if (binr>=0 && binr < pdata->rbins && binphi>=0 && binphi < pdata->phibins && binz>=0 && binz < pdata->zbins) {
      bin_volume=PI*((pdata->minr+(binr+1)*rbinsize)*(pdata->minr+(binr+1)*rbinsize) - (pdata->minr+(binr)*rbinsize)*(pdata->minr+(binr)*rbinsize)) *zbinsize * phibinsize/2/PI;
      A[binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz] += 1./bin_volume;
    } 
  }
  return 0;
}

int observable_radial_flux_density_profile(void* pdata_, double* A, unsigned int n_A) {
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) pdata_;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_r, v_phi, v_z;
    
  for (unsigned i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]/time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]/time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]/time_step;
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    transform_to_cylinder_coordinates(ppos[0]-pdata->center[0], ppos[1]-pdata->center[1], ppos[2]-pdata->center[2], &r, &phi, &z);
    binr  =(int)floor((r-pdata->minr)/rbinsize);
    binphi=(int)floor((phi-pdata->minphi)/phibinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);

    if (binr>=0 && binr < pdata->rbins && binphi>=0 && binphi < pdata->phibins && binz>=0 && binz < pdata->zbins) {
      bin_volume=PI*((pdata->minr+(binr+1)*rbinsize)*(pdata->minr+(binr+1)*rbinsize) - (pdata->minr+(binr)*rbinsize)*(pdata->minr+(binr)*rbinsize)) *zbinsize * phibinsize/2/PI;
      v_r = 1/r*((ppos[0]-pdata->center[0])*v[0] + (ppos[1]-pdata->center[1])*v[1]); 
      v_phi = 1/r/r*((ppos[0]-pdata->center[0])*v[1]-(ppos[1]-pdata->center[1])*v[0]);
      v_z = v[2];
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 0] += v_r/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 1] += v_phi/bin_volume;
      A[3*(binr*pdata->phibins*pdata->zbins + binphi*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}

int observable_flux_density_profile(void* pdata_, double* A, unsigned int n_A) {
  int binx, biny, binz;
  double ppos[3];
  double x, y, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  profile_data* pdata;
  pdata=(profile_data*) pdata_;
  ids=pdata->id_list;
  double xbinsize=(pdata->maxx - pdata->minx)/pdata->xbins;
  double ybinsize=(pdata->maxy - pdata->miny)/pdata->ybins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_x, v_y, v_z;
    
  for (unsigned i = 0; i< n_A; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]*time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]*time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]*time_step;
    memcpy(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memcpy(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    x=ppos[0];
    y=ppos[1];
    z=ppos[2];
    binx  =(int)floor((x-pdata->minx)/xbinsize);
    biny  =(int)floor((y-pdata->miny)/ybinsize);
    binz  =(int)floor((z-pdata->minz)/zbinsize);


    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      bin_volume=xbinsize*ybinsize*zbinsize;
      v_x=v[0];
      v_y=v[1];
      v_z=v[2];
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 0] += v_x/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 1] += v_y/bin_volume;
      A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + 2] += v_z/bin_volume;
    } 
  }
  return 0;
}

int observable_particle_positions(void* idlist, double* A, unsigned int n_A) {
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].r.p[0];
      A[3*i + 1] = partCfg[ids->e[i]].r.p[1];
      A[3*i + 2] = partCfg[ids->e[i]].r.p[2];
  }
  return 0;
}

int observable_particle_forces(void* idlist, double* A, unsigned int n_A) {
  IntList* ids;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  ids=(IntList*) idlist;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].f.f[0]/time_step/time_step*2;
      A[3*i + 1] = partCfg[ids->e[i]].f.f[1]/time_step/time_step*2;
      A[3*i + 2] = partCfg[ids->e[i]].f.f[2]/time_step/time_step*2;
  }
  return 0;
}


int observable_stress_tensor(void* params_p, double* A, unsigned int n_A) {
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  observable_compute_stress_tensor(1,A,n_A);
  return 0;
}


int observable_stress_tensor_acf_obs(void* params_p, double* A, unsigned int n_A) {
  double stress_tensor[9];
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  observable_compute_stress_tensor(1,stress_tensor,9);
  A[0]=stress_tensor[1];
  A[1]=stress_tensor[5];
  A[2]=stress_tensor[6];
  A[3]=stress_tensor[0]-stress_tensor[4];
  A[4]=stress_tensor[0]-stress_tensor[8];
  A[5]=stress_tensor[4]-stress_tensor[8];
  return 0;
}

int observable_structure_factor(void* params_p, double* A, unsigned int n_A) {
  // FIXME Currently scattering length is hardcoded as 1.0
  int l;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr; 
//  DoubleList *scattering_length;
  observable_sf_params* params;
  params = (observable_sf_params*)params_p;
//  scattering_length = params->scattering_length;
  const double scattering_length=1.0;
  order = params->order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }

    for(unsigned int p=0; p<n_A; p++) {
       A[p]   = 0.0;
    }

    l=0;
    //printf("n_A: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
    for (int i=-order; i<=order; i++) {
      for (int j=-order; j<=order; j++) {
        for (int k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
            //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
	    for(int p=0; p<n_part; p++) {
	      qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	      C_sum+= scattering_length * cos(qr);
	      S_sum-= scattering_length * sin(qr);
	    }
            A[l]   =C_sum;
            A[l+1] =S_sum;
            l=l+2;
	  }
	}
      }
    }
    //printf("finished calculating sf\n"); fflush(stdout);
    return 0;
}

int observable_interacts_with (void* params_p, double* A, unsigned int n_A) {
  iw_params *params=(iw_params*)params_p;
  IntList* ids1;
  IntList* ids2;
  int i,j;
//  double dx,dy,dz;
  double dist2;
  double cutoff2=params->cutoff*params->cutoff;
  double pos1[3], pos2[3], dist[3];
  ids1=params->ids1;
  ids2=params->ids2;
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt,"{094 could not sort partCfg} ");
    return -1;
  }
  for ( i = 0; i<ids1->n; i++ ) {
    if (ids1->e[i] >= n_part)
      return 1;
    pos1[0]=partCfg[ids1->e[i]].r.p[0];
    pos1[1]=partCfg[ids1->e[i]].r.p[1];
    pos1[2]=partCfg[ids1->e[i]].r.p[2];
    for ( j = 0; j<ids2->n; j++ ) {
      if (ids2->e[j] >= n_part)
        return 1;
      if (ids2->e[j] == ids1->e[i]) // do not count self-interaction :-)
        continue;
      A[i] = 0;
      pos2[0]=partCfg[ids2->e[j]].r.p[0];
      pos2[1]=partCfg[ids2->e[j]].r.p[1];
      pos2[2]=partCfg[ids2->e[j]].r.p[2];
      get_mi_vector(dist,pos1,pos2);
      dist2= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2];
      if(dist2<cutoff2) {
        A[i] = 1;
	break;
	// interaction found for i, go for next
      }
    }
  }
  return 0;
}





