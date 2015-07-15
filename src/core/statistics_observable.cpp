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
#include "statistics_correlation.hpp"
#include "particle_data.hpp"
#include "integrate.hpp"
#include "lb.hpp"
#include "pressure.hpp"
#include "rotation.hpp"

observable** observables = 0;
int n_observables = 0; 
int observables_autoupdate = 0;

void observable_init(observable* self) {
  self->type = OBSERVABLE;
  self->last_update = 0;
  self->autoupdate = 0;
  self->autoupdate_dt = 0;
}

int observable_calculate(observable* self) {
  int temp = 0;
  if (self->calculate!=0)
    temp=(self->calculate)(self);
  self->last_update = sim_time;
  return temp;
}

int observable_update(observable* self) {
  int temp = 0;
  if (self->update!=0)
    temp=(self->update)(self);
  self->last_update = sim_time;
  return temp;
}

//void write_double(FILE * fp, const double * data, unsigned int n, bool binary)

int observable_write(char *filename, observable *self, bool binary) {
  FILE *f = fopen(filename, "w");
  if(f) {
    unsigned int un;
    /* For stateless observables only the current value is meaningful. */
    if(self->type == OBSERVABLE)
      observable_calculate(self);
    switch(self->type) {
    case AVERAGE:
      write_uint(f, &((observable_average_container*)self->container)->n_sweeps, 1, binary);
    case OBSERVABLE:
      un = self->n;
      write_uint(f, &un, 1, binary);
      write_double(f, self->last_value, self->n, binary);
      fclose(f);
      break;
    default:
      fclose(f);
      return ES_ERROR;
    }
  }
  return ES_ERROR;
}

int observable_read(char *filename, observable *self, bool binary) {
  FILE *f = fopen(filename, "r");
  if(f && !feof(f)) {
    unsigned int un;
    switch(self->type) {
    case AVERAGE:
      read_uint(f, &((observable_average_container*)self->container)->n_sweeps, 1, binary);
    case OBSERVABLE:
      read_uint(f, &un, 1, binary);
      if(self->n != (int)(un))
	return ES_ERROR;
      read_double(f, (double *)self->last_value, self->n, binary);
      fclose(f);
      return ES_OK;
      break;
    default:
      fclose(f);
      return ES_ERROR;
    }
  }
  return ES_ERROR;
}

int observable_calc_particle_velocities(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    A[3*i + 0] = partCfg[ids->e[i]].m.v[0]/time_step;
    A[3*i + 1] = partCfg[ids->e[i]].m.v[1]/time_step;
    A[3*i + 2] = partCfg[ids->e[i]].m.v[2]/time_step;
  }
  return 0;
}

int observable_calc_particle_body_velocities(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;

#ifdef ROTATION

    double RMat[9];
    double vel_lab[3];
    double vel_body[3];

    vel_lab[0] = partCfg[ids->e[i]].m.v[0]/time_step;
    vel_lab[1] = partCfg[ids->e[i]].m.v[1]/time_step;
    vel_lab[2] = partCfg[ids->e[i]].m.v[2]/time_step;
    define_rotation_matrix(&partCfg[ids->e[i]], RMat);

    vel_body[0] = RMat[0 + 3*0]*vel_lab[0] + RMat[0 + 3*1]*vel_lab[1] + RMat[0 + 3*2]*vel_lab[2];
    vel_body[1] = RMat[1 + 3*0]*vel_lab[0] + RMat[1 + 3*1]*vel_lab[1] + RMat[1 + 3*2]*vel_lab[2];
    vel_body[2] = RMat[2 + 3*0]*vel_lab[0] + RMat[2 + 3*1]*vel_lab[1] + RMat[2 + 3*2]*vel_lab[2];

    A[3*i + 0] = vel_body[0];
    A[3*i + 1] = vel_body[1];
    A[3*i + 2] = vel_body[2];

#else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

#endif

  }
  return 0;
}

int observable_calc_particle_angular_momentum(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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

int observable_calc_particle_body_angular_momentum(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  for ( int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;

#ifdef ROTATION

    A[3*i + 0] = partCfg[ids->e[i]].m.omega[0];
    A[3*i + 1] = partCfg[ids->e[i]].m.omega[1];
    A[3*i + 2] = partCfg[ids->e[i]].m.omega[2];

#else

    A[3*i + 0] = 0.0;
    A[3*i + 1] = 0.0;
    A[3*i + 2] = 0.0;

#endif

  }
  return 0;
}

#ifdef ELECTROSTATICS
int observable_calc_particle_currents(observable* self) {
  double* A = self->last_value;
  double charge;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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

int observable_calc_currents(observable* self) {
  double* A = self->last_value;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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

int observable_calc_dipole_moment(observable* self) {
  double* A = self->last_value;
  double charge;
  double j[3] = {0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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

int observable_calc_com_velocity(observable* self) {
  double* A = self->last_value;
  double v_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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

int observable_calc_blocked_com_velocity(observable* self) {
  double* A = self->last_value;
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  n_blocks=self->n/3; 
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

int observable_calc_blocked_com_position(observable* self) {
  double* A = self->last_value;
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  n_blocks=self->n/3; 
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

int observable_calc_com_position(observable* self) {
  double* A = self->last_value;
  double p_com[3] = { 0. , 0., 0. } ;
  double total_mass = 0;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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


int observable_calc_com_force(observable* self) {
  double* A = self->last_value;
  double f_com[3] = { 0. , 0., 0. } ;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
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


int observable_calc_blocked_com_force(observable* self) {
  double* A = self->last_value;
  unsigned int i;
  unsigned int block;
  unsigned int n_blocks;
  unsigned int blocksize;
  unsigned int id;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  n_blocks=self->n/3; 
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


int observable_calc_density_profile(observable* self) {
  double* A = self->last_value;
  int binx, biny, binz;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;

  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  pdata=(profile_data*) self->container;
  ids=pdata->id_list;
  double bin_volume=(pdata->maxx-pdata->minx)*(pdata->maxy-pdata->miny)*(pdata->maxz-pdata->minz)/pdata->xbins/pdata->ybins/pdata->zbins;
    
  for ( int i = 0; i<self->n; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
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

int observable_calc_force_density_profile(observable* self) {
  double* A = self->last_value;
  int binx, biny, binz;
  double ppos[3];
  int img[3];
  IntList* ids;
  profile_data* pdata;

  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  pdata=(profile_data*) self->container;
  ids=pdata->id_list;
  double bin_volume=(pdata->maxx-pdata->minx)*(pdata->maxy-pdata->miny)*(pdata->maxz-pdata->minz)/pdata->xbins/pdata->ybins/pdata->zbins;
    
  for ( int i = 0; i<self->n; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    binx= (int) floor( pdata->xbins*  (ppos[0]-pdata->minx)/(pdata->maxx-pdata->minx));
    biny= (int) floor( pdata->ybins*  (ppos[1]-pdata->miny)/(pdata->maxy-pdata->miny));
    binz= (int) floor( pdata->zbins*  (ppos[2]-pdata->minz)/(pdata->maxz-pdata->minz));
    if (binx>=0 && binx < pdata->xbins && biny>=0 && biny < pdata->ybins && binz>=0 && binz < pdata->zbins) {
      for(int dim = 0; dim < 3; dim++)
        A[3*(binx*pdata->ybins*pdata->zbins + biny*pdata->zbins + binz) + dim] += partCfg[ids->e[i]].f.f[dim]/bin_volume;
    } 
  }
  return 0;
}

#ifdef LB
int observable_calc_lb_velocity_profile(observable* self) {
  double* A= self->last_value;
  void* pdata_ = self->container;
  unsigned int n_A = self->n;
  unsigned int maxi, maxj, maxk;
  double xoffset, yoffset, zoffset;
  double x_incr, y_incr, z_incr;
  double p[3], v[3];
  int linear_index;
  profile_data* pdata;
  pdata=(profile_data*) self->container;


#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU)
    return statistics_observable_lbgpu_velocity_profile((profile_data*) pdata_, A, n_A);
#endif
  if (lattice_switch & LATTICE_LB) {
    for ( int i = 0; i<self->n; i++ ) {
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
    unsigned int i, j, k;
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
  
    for ( int i = 0; i<self->n; i++ ) {
      A[i]*=normalization_factor;
    }

  
  }
  return 0;
}
#endif

#ifdef LB
int observable_calc_lb_radial_velocity_profile(observable* self) {
  double* A = self->last_value;
  void* pdata = self->container;
  unsigned int n_A = self->n;

#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU)
    return statistics_observable_lbgpu_radial_velocity_profile((radial_profile_data*) pdata, A, n_A);
#endif
  
  if (!(lattice_switch & LATTICE_LB))
    return ES_ERROR;

  if (n_nodes==1) {
    mpi_observable_lb_radial_velocity_profile_parallel(pdata, A, n_A);
    return ES_OK;
  } else {
    mpi_observable_lb_radial_velocity_profile();
    MPI_Bcast(pdata, sizeof(radial_profile_data), MPI_BYTE, 0, comm_cart);
    double* data = (double*) malloc(n_A*sizeof(double));
    mpi_observable_lb_radial_velocity_profile_parallel(pdata, data, n_A);
    MPI_Reduce(data, A, n_A, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
    free(data);
    return ES_OK;
  }
}

void mpi_observable_lb_radial_velocity_profile_slave_implementation() {
  radial_profile_data pdata;
  MPI_Bcast(&pdata, sizeof(radial_profile_data), MPI_BYTE, 0, comm_cart);
  unsigned int n_A=3*pdata.rbins*pdata.phibins*pdata.zbins;
  double* data = (double*) malloc(n_A*sizeof(double));
  mpi_observable_lb_radial_velocity_profile_parallel(&pdata, data, n_A);
  MPI_Reduce(data, 0, n_A, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
  free(data);

}

int mpi_observable_lb_radial_velocity_profile_parallel(void* pdata_, double* A, unsigned int n_A) {
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
        if (   p[0] < my_left[0] || p[0]>my_right[0] 
            || p[1] < my_left[1] || p[1]>my_right[1] 
            || p[2] < my_left[2] || p[2]>my_right[2] )
          continue;

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

int observable_calc_radial_density_profile(observable* self) {
  double* A = self->last_value;
  int binr, binphi, binz;
  double ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;
  
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) self->container;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
    
  for (int i = 0; i< self->n; i++ ) {
    A[i]=0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
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

int observable_calc_radial_flux_density_profile(observable* self) {
  double* A = self->last_value;
  int binr, binphi, binz;
  double ppos[3];
  double unfolded_ppos[3];
  double r, phi, z;
  int img[3];
  double bin_volume;
  IntList* ids;

  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  radial_profile_data* pdata;
  pdata=(radial_profile_data*) self->container;
  ids=pdata->id_list;
  double rbinsize=(pdata->maxr - pdata->minr)/pdata->rbins;
  double phibinsize=(pdata->maxphi - pdata->minphi)/pdata->phibins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_r, v_phi, v_z;

  if (self->last_update==sim_time) {
    return ES_ERROR;
  }
    
  for (int i = 0; i< self->n; i++ ) {
    A[i]=0;
  }
  double* old_positions=(double*) pdata->container;
  if (old_positions[0] == CONST_UNITITIALIZED) {
    for (int i = 0; i<ids->n; i++ ) {
      memmove(unfolded_ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
      memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
      unfold_position(unfolded_ppos, img);
      old_positions[3*i+0]=unfolded_ppos[0];
      old_positions[3*i+1]=unfolded_ppos[1];
      old_positions[3*i+2]=unfolded_ppos[2];
    }
    return 0;
  }
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
/* We use folded coordinates here */
    memmove(unfolded_ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    unfold_position(unfolded_ppos, img);
    v[0]=(unfolded_ppos[0] - old_positions[3*i+0]);
    v[1]=(unfolded_ppos[1] - old_positions[3*i+1]);
    v[2]=(unfolded_ppos[2] - old_positions[3*i+2]);
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    // The position of the particle is by definition the middle of old and new position
    ppos[0]+=0.5*v[0]; ppos[1]+=0.5*v[1]; ppos[2]+=0.5*v[2];
    fold_position(ppos, img);
    v[0]/=(sim_time - self->last_update);
    v[1]/=(sim_time - self->last_update);
    v[2]/=(sim_time - self->last_update);
    if (i==0) {
//      printf("(%3.4f) %f %f %f\n", sim_time-self->last_update, v[2], partCfg[ids->e[i]].m.v[2]/time_step,v[2]* partCfg[ids->e[i]].m.v[2]/time_step/time_step);
//      printf("(%3.3f) %f %f", sim_time, old_positions[3*i+2], unfolded_ppos[2]);
    }
    old_positions[3*i+0]=unfolded_ppos[0];
    old_positions[3*i+1]=unfolded_ppos[1];
    old_positions[3*i+2]=unfolded_ppos[2];
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

int observable_calc_flux_density_profile(observable* self) {
  double* A = self->last_value;
  int binx, biny, binz;
  double ppos[3];
  double x, y, z;
  int img[3];
  double bin_volume;
  IntList* ids;

  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  profile_data* pdata;
  pdata=(profile_data*) self->container;
  ids=pdata->id_list;
  double xbinsize=(pdata->maxx - pdata->minx)/pdata->xbins;
  double ybinsize=(pdata->maxy - pdata->miny)/pdata->ybins;
  double zbinsize=(pdata->maxz - pdata->minz)/pdata->zbins;
  double v[3];
  double v_x, v_y, v_z;
    
  for (int i = 0; i< self->n; i++ ) {
    A[i]=0;
  }

  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
    /* We use folded coordinates here */
    v[0]=partCfg[ids->e[i]].m.v[0]/time_step;
    v[1]=partCfg[ids->e[i]].m.v[1]/time_step;
    v[2]=partCfg[ids->e[i]].m.v[2]/time_step;
    memmove(ppos, partCfg[ids->e[i]].r.p, 3*sizeof(double));
    memmove(img, partCfg[ids->e[i]].l.i, 3*sizeof(int));
    fold_position(ppos, img);
    // The position of the particle is by definition the middle of old and new position
  
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

int observable_calc_particle_positions(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].r.p[0];
      A[3*i + 1] = partCfg[ids->e[i]].r.p[1];
      A[3*i + 2] = partCfg[ids->e[i]].r.p[2];
  }
  return 0;
}

int observable_calc_particle_forces(observable* self) {
  double* A = self->last_value;
  IntList* ids;
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  ids=(IntList*) self->container;
  for (int i = 0; i<ids->n; i++ ) {
    if (ids->e[i] >= n_part)
      return 1;
      A[3*i + 0] = partCfg[ids->e[i]].f.f[0]/time_step/time_step*2;
      A[3*i + 1] = partCfg[ids->e[i]].f.f[1]/time_step/time_step*2;
      A[3*i + 2] = partCfg[ids->e[i]].f.f[2]/time_step/time_step*2;
  }
  return 0;
}


int observable_stress_tensor(observable* self) {
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
    return -1;
  }
  observable_compute_stress_tensor(1,self->last_value,self->n);
  return 0;
}


int observable_calc_stress_tensor_acf_obs(observable* self) {
  double* A = self->last_value;
  double stress_tensor[9];
  if (!sortPartCfg()) {
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
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

int observable_update_average(observable* self) {
    observable_average_container* data = (observable_average_container*) self->container;
    data->n_sweeps++;
    int error = observable_calculate(data->reference_observable);
    if ( error != 0)
      return 1;
    double factor = 1 / (double) data->n_sweeps;
    for (int i =0; i<self->n; i++) {
      self->last_value[i] = (1-factor)*self->last_value[i] + factor*data->reference_observable->last_value[i];
    }
    return 0;
}

int observable_reset_average(observable* self) {
    observable_average_container* data = (observable_average_container*) self->container;
    data->n_sweeps=0;
    int error = observable_calculate(data->reference_observable);
    for (int i =0; i<self->n; i++) {
      self->last_value[i] = 0;
    }
    return error;
}


int observable_calc_structure_factor(observable* self) {
  double* A = self->last_value;
  // FIXME Currently scattering length is hardcoded as 1.0
  int l;
  int order, order2, n;
  double twoPI_L, C_sum, S_sum, qr; 
  //  DoubleList *scattering_length;
  observable_sf_params * params;
  params = (observable_sf_params*) self->container;
  //  scattering_length = params->scattering_length;
  const double scattering_length=1.0;
  order = params->order;
  order2=order*order;
  twoPI_L = 2*PI/box_l[0];
  
  if (!sortPartCfg()) {
    ostringstream msg;
    msg <<"could not sort partCfg";
    runtimeError(msg);
    return -1;
  }

  
  l=0;
  float partCache[n_part*3];
  for(int p=0; p<n_part; p++) {
    for (int i=0;i<3;i++){
      partCache[3*p+i]=partCfg[p].r.p[i];
    }
  }
  //printf("self->n: %d, dim_sf: %d\n",n_A, params.dim_sf); fflush(stdout);
  for(int i=-order; i<=order; i++) {
    for(int j=-order; j<=order; j++) {
      for(int k=-order; k<=order; k++) {
	n = i*i + j*j + k*k;
	if ((n<=order2) && (n>=1)) {
	  C_sum = S_sum = 0.0;
	  //printf("l: %d, n: %d %d %d\n",l,i,j,k); fflush(stdout);
	  for(int p=0; p<n_part; p++) {
	    //qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
	    qr = twoPI_L * ( i*partCache[3*p+0] + j*partCache[3*p+1] + k*partCache[3*p+2] );
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
  l = 0;
  for(int k=0;k<self->n;k++) {
    //devide by the sqrt(number_of_particle) due to complex product and no k-vector averaging so far
    A[k] /= sqrt(n_part);
  }
  //printf("finished calculating sf\n"); fflush(stdout);
  return 0;
}


int observable_calc_structure_factor_fast(observable* self) {
  //printf("calculating\n");
  double* A = self->last_value;
  // FIXME Currently scattering length is hardcoded as 1.0
  observable_sf_params * params = (observable_sf_params*) self->container;
  const int k_max = params->order * params->k_density;
  const double scattering_length=1.0;
  const double twoPI_L = 2*PI/box_l[0];
  
  if (!sortPartCfg()) {
    ostringstream msg;
    msg <<"could not sort partCfg";
    runtimeError(msg);
    return -1;
  }
  
  for(int p=0; p<self->n; p++) {
    A[p]   = 0.0;
  }
  
  float partCache[n_part*3];
  for(int p=0; p<n_part; p++) {
    for (int i=0;i<3;i++){
      partCache[3*p+i]=partCfg[p].r.p[i];
    }
  }
  int k_density = params->k_density;
  int l=0;
  for(int k=0; k<k_max; k++) {
    int order=k/k_density+1;
    switch (k % k_density){
    case 0: // length sqrt(1)
      for (int dir=0;dir<3;dir++){
	double C_sum = 0;
	double S_sum = 0;
	for(int p=0; p<n_part; p++) {
	//double qr = twoPI_L * k  * ( ix*partCache[3*p+0] + iy*partCache[3*p+1] + iz*partCache[3*p+2] );
	  double qr = twoPI_L * order * ( partCache[3*p+dir]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 1: // length sqrt(2)
      for (int dir=0;dir<6;dir++){
	int fac1,fac2,off1,off2;
	switch (dir){
	case 0: fac1= 1; off1=0; fac2= 1; off2=1; break;
	case 1: fac1= 1; off1=0; fac2= 1; off2=2; break;
	case 2: fac1= 1; off1=1; fac2= 1; off2=2; break;
	case 3: fac1=-1; off1=0; fac2= 1; off2=1; break;
	case 4: fac1=-1; off1=0; fac2= 1; off2=2; break;
	case 5: fac1=-1; off1=1; fac2= 1; off2=2; break;
	}
	double C_sum = 0;
	double S_sum = 0;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+off1]*fac1+ partCache[3*p+off2]*fac2);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 2: // length sqrt(3)
      for (int dir=0;dir<4;dir++){
	double C_sum = 0;
	double S_sum = 0;
	int fac1=(1-2*(dir%2));
	int fac2=(1-2*(dir/2));
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+0]*fac1 + partCache[3*p+1]*fac2 + partCache[3*p+2]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 3: // length sqrt(5) 
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3));
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)]*fac1 + partCache[3*p+((dir+1)%3)]*2);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 4: // length sqrt(6) 
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3))*2;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)] + partCache[3*p+((dir+1)%3)]*fac1 + partCache[3*p+((dir+2)%3)]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    case 5: // length sqrt(9)
      for (int dir=0;dir<6;dir++){
	double C_sum = 0;
	double S_sum = 0;
	// pick 6 random vectors
	int fac1=(1-2*(dir/3))*2;
	for(int p=0; p<n_part; p++) {
	  double qr = twoPI_L * order * ( partCache[3*p+(dir%3)]*fac1 + partCache[3*p+((dir+1)%3)]*2 + partCache[3*p+((dir+2)%3)]);
	  C_sum+= scattering_length * cos(qr);
	  S_sum+= scattering_length * sin(qr);
	}
	A[l]   =C_sum;
	A[l+1] =S_sum;
	l+=2;
      }
      break;
    default:
      ostringstream msg;
      msg <<"so many samples per order not yet implemented";
      runtimeError(msg);
      return -1;
    }
  }
  for(int l=0;l<self->n;l++) {
    //devide by the sqrt of number_of_particle, average later
    A[l] /= sqrt(n_part);
  }
  return 0;
}

int observable_calc_interacts_with (observable* self) {
  double* A = self->last_value;
  iw_params *params=(iw_params*) self->container;
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
      ostringstream msg;
      msg <<"could not sort partCfg";
      runtimeError(msg);
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

int observable_calc_rdf(observable* self){
  if (!sortPartCfg()) {
    return 1;
  }
  double * last = self->last_value;
  rdf_profile_data * rdf_data = (rdf_profile_data *) self->container;
  calc_rdf(rdf_data->p1_types, rdf_data->n_p1,
	   rdf_data->p2_types, rdf_data->n_p2,
	   rdf_data->r_min,    rdf_data->r_max,
	   rdf_data->r_bins,   last);
  return 0;
}


int observable_radial_density_distribution(observable* self){
  if (!sortPartCfg()) {
    ostringstream errtxt; // = runtime_error(128);
    errtxt << "{094 could not sort partCfg} ";
	runtimeError(errtxt);
    return -1;
  }

  radial_density_data *r_data = (radial_density_data *) self->container;
  IntList *ids;  
  if ( GC_init && Type_array_init ) {
	  ids = (IntList *) malloc(sizeof(IntList));

	  //using the grandcanonical scheme, always update the particle id list
	  ids->e = (int *) malloc(sizeof(int)*type_array[Index.type[r_data->type]].max_entry);
	  memmove(ids->e, type_array[Index.type[r_data->type]].id_list, type_array[Index.type[r_data->type]].max_entry*sizeof(int));
	  ids->n = type_array[Index.type[r_data->type]].max_entry;
	  ids->max = type_array[Index.type[r_data->type]].cur_size;
  } else { 
	  ids = r_data->id_list;
  }
  //r_data->id_list = ids;
  //ids = r_data->id_list;

  double* A = self->last_value;
  int n_A   = self->n;
  double start_point[3];
  double end_point[3];
  int image_box[3];
  if ( r_data->id_flag ) {
	  // Using particle_ids to specify the start and endpoints
	  memmove(start_point, partCfg[r_data->start_point_id].r.p, 3*sizeof(double));
	  memmove(image_box, partCfg[r_data->start_point_id].l.i, 3*sizeof(int));
	  unfold_position(start_point, image_box);
	  memmove(end_point, partCfg[r_data->end_point_id].r.p, 3*sizeof(double));
	  memmove(image_box, partCfg[r_data->end_point_id].l.i, 3*sizeof(int));
	  unfold_position(end_point, image_box);
  } else {
	  memmove(start_point, r_data->start_point, 3*sizeof(double));
	  memmove(end_point, r_data->end_point, 3*sizeof(double));
  }

  double *bin_volume = (double *) malloc(sizeof(double)*r_data->rbins);
 
  double part_pos[3];
  double AB[3];		// normalized normal vector pointing to start point
  double BA[3];		// ...									end point
  double Apart[3];	// vector difference start_point - part_pos
  double norm_Apart;
  double normAB;
  double dist;
  double angle;
  double r_dist;
  double tmp[3];
  double tmp2[3];
  char periodic_points = 0;
  int frac = 0;
  double epsilon = 1e-9;

  double rbin_size = (r_data->maxr - r_data->minr)/(r_data->rbins);
  int rbin_id;

  if ( r_data->id_flag ) {
	  // in the case of particles use minimal image convention for axis definition
	  get_mi_vector(AB, start_point, end_point);
	  get_mi_vector(BA, end_point, start_point);
	  if ( normr(AB) < epsilon )  {
		  //start and end point are periodic images of one another
		  periodic_points = 1;
	  }
	  end_point[0] = start_point[0] - AB[0];
	  end_point[1] = start_point[1] - AB[1];
	  end_point[2] = start_point[2] - AB[2];
  } else {
	  // otherwise use the given absolute positions but check whether those points are too close
	  get_mi_vector(AB, start_point, end_point);
	  if ( normr(AB) < epsilon )  {
		  //start and end point are periodic images of one another
		  periodic_points = 1;
	  }
	  vecsub(start_point, end_point, AB);
	  vecsub(end_point, start_point, BA);
  }

  //printf("id: %d, lenAB: %1.3e, lenBA: %1.3e, %f %f %f, %f %f %f\n", r_data->id_flag, normr(AB), normr(BA), start_point[0], start_point[1], start_point[2], end_point[0], end_point[1], end_point[2]);
  normAB = normr(AB); 
  for (int i=0; i<r_data->rbins; i++) {
	  bin_volume[i] = normAB*( pow(r_data->minr + (i+1)*rbin_size, 2) - pow(r_data->minr + i*rbin_size, 2) ) * PI;
  }

  unit_vector(AB, AB);
  unit_vector(BA, BA);

  for (int i=0; i<n_A; i++)
	  A[i] = 0.0;

  for (int i=0; i < ids->n; i++){
	  part_pos[0]=partCfg[ids->e[i]].r.p[0];
	  part_pos[1]=partCfg[ids->e[i]].r.p[1];
	  part_pos[2]=partCfg[ids->e[i]].r.p[2];
	  // that might seem weird, but this ensures that the used particle position is the
	  // closest image to one of the two points that define the axis
	  get_mi_vector(tmp, part_pos, start_point);
	  get_mi_vector(tmp2, part_pos, end_point);

	  if ( periodic_points ) { 
		  // the axis spans the whole box, so any particle is in 
		  // scope
		  dist = -1;
	  } else if ( normr(tmp) < normr(tmp2) ) {
		  part_pos[0] = start_point[0] + tmp[0];
		  part_pos[1] = start_point[1] + tmp[1];
		  part_pos[2] = start_point[2] + tmp[2];
		  dist = scalar(AB, tmp);
	  } else {
		  part_pos[0] = end_point[0] + tmp2[0];
		  part_pos[1] = end_point[1] + tmp2[1];
		  part_pos[2] = end_point[2] + tmp2[2];
		  dist =  scalar(BA, tmp2);
	  }

	  if (dist > 0 ) { 

		  continue;
	  }

	  // particle in scope
	  // now calculate the distance from the given axis using simple geometry
	  vecsub(start_point, part_pos, Apart);

	  norm_Apart = normr(Apart);
	  angle = acos(scalar(Apart, AB)/norm_Apart);
	  r_dist= sin(angle)*norm_Apart;

	  rbin_id = (int) floor( (r_dist - r_data->minr)/rbin_size );
	  if ( rbin_id >= r_data->rbins || rbin_id < 0 )
		  continue;

	  A[rbin_id]+=1.0/bin_volume[rbin_id];
	  frac++;
  }
//  printf("fraction of parts: %d %d\n", frac, ids->n);
  free(bin_volume);
  if ( GC_init && Type_array_init ) {
	  free(ids->e);
	  free(ids);
	}
  return 0;
}

int observable_spatial_polymer_properties(observable* self){
    if (!sortPartCfg()) {
      ostringstream errtxt ; //= runtime_error(128);
	  errtxt << "{094 could not sort partCfg}";
	  runtimeError(errtxt);
      return -1;
    }
	double* A = self->last_value;
	int poly_len = self->n;
	for (int i = 0; i<poly_len; i++ )
		A[i]=0.0;

#ifdef ELECTROSTATICS
	spatial_polym_data *p_data = (spatial_polym_data *) self->container;
	IntList *ids=p_data->id_list;
	for (int i = 0; i<ids->n; i++){
		A[i%poly_len] += partCfg[ids->e[i]].p.q/( (double) p_data->npoly);
	}
#endif
	return 0;
}

int observable_persistence_length(observable* self){
    if (!sortPartCfg()) {
      ostringstream errtxt; // = runtimeError(128);
      errtxt << "{094 could not sort partCfg}";
	  runtimeError(errtxt);
      return -1;
    }
	double* A = self->last_value;
	spatial_polym_data *p_data = (spatial_polym_data *) self->container;
	IntList *ids=p_data->id_list;

	double v1[3];
	double v2[3];
	double abs1, abs2;
	int num_parts = ids->n;
	int cut_off = p_data->cut_off;
	int n_A = self->n;

	for (int i = 0; i<n_A; i++ )
		A[i]=0.0;

	for (int i = 0; i<n_A; i++) {
		// i : distance between polymer segments
		for ( int j=cut_off; j<num_parts - i -cut_off; j++ ) {
			vecsub(partCfg[ids->e[j]].r.p, partCfg[ids->e[j+1]].r.p, v1);
			abs1 = normr(v1);
			vecsub(partCfg[ids->e[i+j]].r.p, partCfg[ids->e[i+j+1]].r.p, v2);
			abs2 = normr(v2);
			A[i] += scalar(v1, v2)/(abs1*abs2 * (double) (num_parts - i - 2*cut_off) );
		}
	}
	return 0;
}	

	

int observable_polymer_k_distribution(observable *self){
	if (!sortPartCfg()) {
		ostringstream errtxt;
		errtxt << "{094 could not sort partCfg}";
		runtimeError(errtxt);
		return -1;
	}
	double* A = self->last_value;
	for (int i = 0; i < self->n; i++) {
		A[i] = 0.0;
	}
	k_dist_data *data = (k_dist_data *) self->container;
	IntList *ids = data->id_list;
	int npoly = data->npoly;
	int poly_len = data-> poly_len;
	int k = data->k;
	double dist_vec[3];
	double dist;
	double r_min = data->r_min;
	double r_max = data->r_max;
	int n_bins   = data->n_bins;
	int bin_id =0;
	double bin_size = (r_max - r_min) / n_bins;
	int number_of_pairs =(int) ( floor(poly_len / (k +1 ) ) + ( (poly_len - 1)%(k+1) == 0 ? 1 : 0 ) - 1);
	for (int i = 0; i < npoly; i++) 
		for (int j = 0; j < poly_len - k; j+= k) {
			get_mi_vector(dist_vec, partCfg[ids->e[i*poly_len + j]].r.p, partCfg[ids->e[i*poly_len + j + k]].r.p);
			dist = normr(dist_vec);
			bin_id = (int) floor( (dist - r_min)/bin_size );
			if (bin_id < n_bins && bin_id >= 0) {
				A[bin_id] += 1.0/(npoly * number_of_pairs);
			}
		}
	return 0;
}	

void autoupdate_observables() {
  int i;
  for (i=0; i<n_observables; i++) {
//    printf("checking observable %d autoupdate is %d \n", i, observables[i]->autoupdate);
    if (observables[i]->autoupdate && sim_time-observables[i]->last_update>observables[i]->autoupdate_dt*0.99999) {
//      printf("updating %d\n", i);
      observable_update(observables[i]);
    }
  }
}
