#include "particle_data.hpp"
#include "utils.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "cells.hpp"
#include "initialize.hpp" 
#include "integrate.hpp" 


void local_rotate_system(double phi, double theta, double alpha)
{
  // Culculate center of mass
  double com[3];
  int N=0; // Num of particles

  int c, np, i;
  Particle *part;
  Cell *cell;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      for (int j=0;j<3;j++){
        com[j]+=part[i].r.p[j];
      }
      N++;
    }
  }
  for (int j=0;j<3;j++) com[j]/=N;

  
  // Rotation axis in carthesian coordinates
  double axis[3];
  axis[0]=sin(theta)*cos(phi);
  axis[1]=sin(theta)*sin(phi);
  axis[2]=cos(theta);

  // Rotate particle coordinates
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np = cell->n;
        
    for(i = 0; i < np; i++) {
      // Move the center of mass of the system to the origin
      for (int j=0;j<3;j++){
        part[i].r.p[j]-=com[j];
      }
      // Rotate
      double res[3];
      vec_rotate(axis,alpha,part[i].r.p,res);
      // Write back result and shift back the center of mass
      for (int j=0;j<3;j++){
        part[i].r.p[j]=com[j]+res[j];
      }
    }
  }


  resort_particles =1;
  announce_resort_particles();

}


/** Rotate all particle coordinates around an axis given by phi,theta through the center of mass by an angle alpha */
void rotate_system(double phi, double theta, double alpha)
{
  if (n_nodes!=1) {
      ostringstream msg;
      msg <<"Rotate_system only works on a single cpu core";
      runtimeError(msg);
  }

  local_rotate_system(phi, theta, alpha);
}

