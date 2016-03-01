/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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
#include "external_potential.hpp"
#include "lattice.hpp"
#include "communication.hpp"
#include "integrate.hpp"

ExternalPotential* external_potentials;
int n_external_potentials;

void external_potential_pre_init() {
  external_potentials = NULL;
  n_external_potentials = 0;
}


int generate_external_potential(ExternalPotential** e) {
  external_potentials = (ExternalPotential*) Utils::realloc(external_potentials,
		  (n_external_potentials+1) * sizeof(ExternalPotential));
  *e = &external_potentials[n_external_potentials];
  n_external_potentials++;
  (*e)->energy = 0;


//  e = &external_potentials[n_external_potentials-1];
  return ES_OK;
}     

int external_potential_tabulated_init(int number, char* filename, int n_particle_types, double* scale) {
  ExternalPotentialTabulated* e = &external_potentials[number].tabulated;

  if (strlen(filename)>MAX_FILENAME_SIZE)
    return ES_ERROR;
  strcpy((char*)&(e->filename), filename);
  external_potentials[number].type=EXTERNAL_POTENTIAL_TYPE_TABULATED;
  external_potentials[number].scale = (double*) Utils::malloc(n_particle_types*sizeof(double));
  external_potentials[number].n_particle_types = n_particle_types;
  for (int i = 0; i < n_particle_types; i++) {
    external_potentials[number].scale[i]=scale[i];
  }
  mpi_external_potential_broadcast(number);
  mpi_external_potential_tabulated_read_potential_file(number);
  return ES_OK;
}

int lattice_read_file(Lattice* lattice, char* filename);

int external_potential_tabulated_read_potential_file(int number) {
  return lattice_read_file(&(external_potentials[number].tabulated.potential),
		  external_potentials[number].tabulated.filename);
}

int lattice_read_file(Lattice* lattice, char* filename) {
 // ExternalPotentialTabulated *e = &(external_potentials[number].e.tabulated);
  FILE* infile = fopen(filename, "r");
  
  if (!infile)  {
    runtimeErrorMsg() <<"Could not open file "<< filename << "\n";
    return ES_ERROR;
  }
  char first_line[100];
  char* token;
  double res[3];
  double size[3];
  double offset[3]={0,0,0};
  int dim=0;
  if (fgets(first_line, 100, infile) == NULL) {
      fprintf(stderr, "Nothing read from file\n");
      return ES_ERROR;
  }

  token = strtok(first_line, " \t");
  if (!token) { fprintf(stderr, "Error reading dimensionality\n"); return ES_ERROR; }
  dim = atoi(token);
  if (dim<=0)  { fprintf(stderr, "Error reading dimensionality\n"); return ES_ERROR; }
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[0]\n"); return ES_ERROR; }
  size[0] = atof(token);

  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[1]\n"); return ES_ERROR; }
  size[1] = atof(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[2]\n"); return ES_ERROR;}
  size[2] = atof(token);

  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read res[0]\n"); return ES_ERROR;}
  res[0] = atof(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read res[1]\n"); return ES_ERROR;}
  res[1] = atof(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read res[2]\n"); return ES_ERROR;}
  res[2] = atof(token);

  token = strtok(NULL, " \t");
  if (token) {
    offset[0]=atof(token);
    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read offset[1]\n"); return ES_ERROR;}
    offset[1] = atof(token);
    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read offset[2]\n"); return ES_ERROR;}
    offset[2] = atof(token);
  }
  lattice->offset[0]=offset[0];
  lattice->offset[1]=offset[1];
  lattice->offset[2]=offset[2];

  int halosize=1;

  if (size[0] > 0 && fabs(size[0] - box_l[0]) > ROUND_ERROR_PREC) {
    runtimeErrorMsg() <<"Box size in x is wrong "<< size[0] << " vs " << box_l[0] <<"\n";
    return ES_ERROR;
  }
  if (size[1] > 0 && fabs(size[1] - box_l[1]) > ROUND_ERROR_PREC) {
    runtimeErrorMsg() <<"Box size in y is wrong "<< size[1] << " vs " << box_l[1] <<"\n";
    return ES_ERROR;
  }
  if (size[2] > 0 && fabs(size[2] - box_l[2]) > ROUND_ERROR_PREC) {
    runtimeErrorMsg() <<"Box size in z is wrong "<< size[2] << " vs " << box_l[2] <<"\n";
    return ES_ERROR;
  }


  if (res[0] > 0)
    if (skin/res[0]>halosize) halosize = (int)ceil(skin/res[0]);
  if (res[1] > 0)
    if (skin/res[1]>halosize) halosize = (int)ceil(skin/res[1]);
  if (res[2] > 0)
    if (skin/res[2]>halosize) halosize = (int)ceil(skin/res[2]);

  // Now we count how many entries we have:

  lattice->init(res, offset, halosize, dim);
  lattice->interpolation_type = INTERPOLATION_LINEAR;

  char* line = (char*) Utils::malloc((3+dim)*ES_DOUBLE_SPACE);
  double pos[3];
  double f[3];
  int i;
  
  while (fgets(line, 200, infile)) {
    if (strlen(line)<2)
      continue;
    token = strtok(line, " \t");
    if (!token) { fprintf(stderr, "Could not read pos[0]\n"); return ES_ERROR; }
    pos[0] = atof(token);

    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read pos[1] in line:\n%s\n", line); return ES_ERROR; }
    pos[1] = atof(token);
    
    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read pos[1]\n"); return ES_ERROR; }
    pos[2] = atof(token);
    for (i=0; i<dim;i++) {
      token = strtok(NULL, " \t");
      if (!token) { fprintf(stderr, "Could not read f[%d]\n", i); return ES_ERROR; }
      f[i] = atof(token);
    }
    lattice->set_data_for_global_position_with_periodic_image(pos, f);
  }
  free(line);

  write_local_lattice_to_file("lattice", lattice);
  
  if (check_runtime_errors()!=0)
    return ES_ERROR;
  return ES_OK;
}


int write_local_lattice_to_file(const char* filename_prefix, Lattice* lattice) {
  index_t index[3];
  double pos[3];
  int i,j,k;
  double *d;

  char filename[60];
  //Lattice* l = lattice;
  sprintf(filename, "%s_%02d.dat", filename_prefix, this_node);
  FILE* outfile = fopen(filename , "w");
  fprintf(outfile,"grid %d %d %d\n", lattice->grid[0], lattice->grid[1], lattice->grid[2]);
  fprintf(outfile,"halo_grid %d %d %d\n", lattice->halo_grid[0], lattice->halo_grid[1], lattice->halo_grid[2]);
  fprintf(outfile,"halo_size %d\n", lattice->halo_size);
  
  fprintf(outfile,"grid_volume %ld\n", lattice->grid_volume);
  fprintf(outfile,"halo_grid_volume %ld\n", lattice->halo_grid_volume);
  fprintf(outfile,"halo_grid_surface %ld\n", lattice->halo_grid_surface);
  fprintf(outfile,"halo_offset %ld\n", lattice->halo_offset);

  fprintf(outfile,"dim %d\n", lattice->dim);

  fprintf(outfile,"agrid %f %f %f\n", lattice->agrid[0], lattice->agrid[1], lattice->agrid[2]);
 
  fprintf(outfile,"offset %f %f %f\n", lattice->offset[0], lattice->offset[1], lattice->offset[2]);
  fprintf(outfile,"local_offset %f %f %f\n", lattice->local_offset[0], lattice->local_offset[1], lattice->local_offset[2]);
  fprintf(outfile,"local_index_offset %d %d %d\n", lattice->local_index_offset[0], lattice->local_index_offset[1], lattice->local_index_offset[2]);


  fprintf(outfile, "element_size %lu\n", lattice->element_size);

  
  for (i=0; i<lattice->halo_grid[0]; i++) 
    for (j=0; j<lattice->halo_grid[1]; j++) 
      for (k=0; k<lattice->halo_grid[2]; k++) {
        index[0]=i; index[1] = j; index[2] = k;
        lattice->get_data_for_halo_index(index, (void**) &d);
        lattice->map_halo_index_to_pos(index, pos);
//        map_local_index_to_pos(&e->lattice, index, pos);
        fprintf(outfile, "%f %f %f %f \n",pos[0], pos[1], pos[2], d[0]);
      } 
  fclose(outfile);
  return ES_OK;
}

void add_external_potential_tabulated_forces(ExternalPotential* e, Particle* p) {
  if (p->p.type >= e->n_particle_types || e->scale[p->p.type] == 0 ) {
    return;
  }
  double field[3];
  double ppos[3];
  int    img[3];
  memmove(ppos, p->r.p, 3*sizeof(double));
  memmove(img, p->r.p, 3*sizeof(int));
  fold_position(ppos, img);
 
  e->tabulated.potential.interpolate_gradient(p->r.p, field);
  p->f.f[0]-=e->scale[p->p.type]*field[0];
  p->f.f[1]-=e->scale[p->p.type]*field[1];
  p->f.f[2]-=e->scale[p->p.type]*field[2];
//  printf("%d %f force: %f %f %f\n", p->p.type, e->scale[p->p.type], e->scale[p->p.type]*field[0], e->scale[p->p.type]*field[1], e->scale[p->p.type]*field[2]);
}

void add_external_potential_forces(Particle* p) {
  for (int i = 0; i < n_external_potentials; i++) {
    if (external_potentials[i].type==EXTERNAL_POTENTIAL_TYPE_TABULATED) {
      add_external_potential_tabulated_forces(&external_potentials[i], p);
    } else {
      runtimeErrorMsg() <<"unknown external potential type";
      return;
    }
  }
}


void add_external_potential_tabulated_energy(ExternalPotential* e, Particle* p) {
  if (p->p.type >= e->n_particle_types) {
    return;
  }
  double potential;
  double ppos[3];
  int img[3];
  memmove(ppos, p->r.p, 3*sizeof(double));
  memmove(img, p->r.p, 3*sizeof(int));
  fold_position(ppos, img);
 
  e->tabulated.potential.interpolate(p->r.p, &potential);
  e->energy += e->scale[p->p.type] * potential;
}

void add_external_potential_energy(Particle* p) {
  for (int i=0; i<n_external_potentials; i++) {
    if (external_potentials[i].type==EXTERNAL_POTENTIAL_TYPE_TABULATED) {
      add_external_potential_tabulated_energy(&external_potentials[i], p);
    } else {
      runtimeErrorMsg() <<"unknown external potential type";
      return;
    }
  }
}

void external_potential_init_energies() {
  for (int i = 0; i<n_external_potentials; i++) {
    external_potentials[i].energy=0;
  }
}

