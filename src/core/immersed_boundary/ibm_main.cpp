/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation

#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "lb.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "halo.hpp"
#include "lb-boundaries.hpp"
#include "immersed_boundary/ibm_main.hpp"
#include "immersed_boundary/ibm_cuda_interface.hpp"

// Dummy functions. They are required by Espresso, but they don't do anything here.
// We have our own update functions.
void update_mol_pos_particle(Particle *) {};
void update_mol_vel_particle(Particle *) {};
void distribute_mol_force() {};

// ****** Functions for internal use ********

void CoupleIBMParticleToFluid(Particle *p);
void ParticleVelocitiesFromLB_CPU();
bool IsHalo(const int indexCheck);
void GetIBMInterpolatedVelocity(double *p, double *const v, double *const forceAdded);

// ***** Internal variables ******

bool *isHaloCache = NULL;

// ******** Variables from other espresso files *****
extern HaloCommunicator update_halo_comm;

/****************
  IBM_ForcesIntoFluid_CPU
 
 Puts the calculated force stored on the ibm particles into the fluid by updating the lbfields structure
 Calls couple_trace_to_fluid for each node
 Called from the integrate loop right after the forces have been calculated
*****************/

void IBM_ForcesIntoFluid_CPU()
{
  
  // Halo has already been sent. Check for safety
  if (lbpar.resend_halo)
  {
    printf("Error. Halo should already be sent!\n");
    exit(1);
  }
  
  // Update the forces on the ghost particles
  ghost_communicator(&cell_structure.ibm_ghost_force_comm);
  
  
  // Loop over local cells
  for (int c = 0; c < local_cells.n; c++)
  {
    Cell *cell = local_cells.cell[c] ;
    Particle *p = cell->part ;
    const int np = cell->n ;
    
    for (int i = 0; i < np; i++)
      if (p[i].p.isVirtual)
        CoupleIBMParticleToFluid(&p[i]);
  }
  
  
  // Loop over ghost cells
  for (int c = 0; c < ghost_cells.n; c++)
  {
    Cell *cell = ghost_cells.cell[c] ;
    Particle *p = cell->part ;
    const int np = cell->n ;
    
    for (int i = 0; i < np; i++)
    {
      // for ghost particles we have to check if they lie
      // in the range of the local lattice nodes
      if (p[i].r.p[0] >= my_left[0]-0.5*lblattice.agrid[0]
          && p[i].r.p[0] < my_right[0]+0.5*lblattice.agrid[0]
          && p[i].r.p[1] >= my_left[1]-0.5*lblattice.agrid[1]
          && p[i].r.p[1] < my_right[1]+0.5*lblattice.agrid[1]
          && p[i].r.p[2] >= my_left[2]-0.5*lblattice.agrid[2]
          && p[i].r.p[2] < my_right[2]+0.5*lblattice.agrid[2])
      {
        
        if (p[i].p.isVirtual)
          CoupleIBMParticleToFluid(&p[i]);
      }
    }
  }
}

/***************
  IBM_ResetLBForces_CPU
Called from the integrate loop directly after the IBM particle update
Usually the reset would be done by Espresso after the LB update. But we need to keep the forces till after the position update for the f/2 term
****************/

void IBM_ResetLBForces_CPU()
{
  for (int i = 0; i<lblattice.halo_grid_volume; ++i)
  {
#ifdef EXTERNAL_FORCES
    // unit conversion: force density
    lbfields[i].force[0] = lbpar.ext_force[0]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau;
    lbfields[i].force[1] = lbpar.ext_force[1]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau;
    lbfields[i].force[2] = lbpar.ext_force[2]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau;
#else
    lbfields[i].force[0] = 0.0;
    lbfields[i].force[1] = 0.0;
    lbfields[i].force[2] = 0.0;
    lbfields[i].has_force = 0;
#endif
  }
}

/*************
  IBM_UpdateParticlePositions
This function is called from the integrate right after the LB update
Interpolates LB velocity at the particle positions and propagates the particles
**************/

void IBM_UpdateParticlePositions()
{
  // Get velocities
  if (lattice_switch & LATTICE_LB) ParticleVelocitiesFromLB_CPU();
#ifdef LB_GPU
  if (lattice_switch & LATTICE_LB_GPU) ParticleVelocitiesFromLB_GPU();
#endif
  
  
  // Do update: Euler
  const double skin2 = SQR(0.5 * skin);
  // Loop over particles in local cells
  for (int c = 0; c < local_cells.n; c++)
  {
    const Cell *const cell = local_cells.cell[c];
    Particle *const p  = cell->part;
    for(int j = 0; j < cell->n; j++)
      if (p[j].p.isVirtual)
      {
        if ( !( p[j].p.ext_flag & 2 ) )
          p[j].r.p[0] = p[j].r.p[0] + p[j].m.v[0]*time_step;
        if ( !( p[j].p.ext_flag & 4 ) )
          p[j].r.p[1] = p[j].r.p[1] + p[j].m.v[1]*time_step;
        if ( !( p[j].p.ext_flag & 8 ) )
          p[j].r.p[2] = p[j].r.p[2] + p[j].m.v[2]*time_step;
        
        // Check if the particle might have crossed a box border (criterion see e-mail Axel 28.8.2014)
        // if possible resort_particles = 1
        const double dist2 = distance2( p[j].r.p, p[j].l.p_old);
        if ( dist2 > skin2 ) { resort_particles = 1; }
      }
  }
  
  // This function spreads the resort_particles variable across the nodes
  // If one node wants to resort, all nodes do it
  announce_resort_particles();
}

/*************
   CoupleIBMParticleToFluid
This function puts the momentum of a given particle into the LB fluid - only for CPU
**************/

void CoupleIBMParticleToFluid(Particle *p)
{
  // Convert units from MD to LB
  double delta_j[3];
  delta_j[0] = p->f.f[0]*time_step*lbpar.tau/lbpar.agrid;
  delta_j[1] = p->f.f[1]*time_step*lbpar.tau/lbpar.agrid;
  delta_j[2] = p->f.f[2]*time_step*lbpar.tau/lbpar.agrid;
  
  // Get indices and weights of affected nodes using discrete delta function
  index_t node_index[8];
  double delta[6];
  lblattice.map_position_to_lattice(p->r.p,node_index,delta);
  
  // Loop over all affected nodes
  for ( int z = 0; z < 2; z++)
  {
    for (int y = 0; y < 2; y++)
    {
      for (int x = 0; x < 2;x++)
      {
        // Do not put force into a halo node
        if ( !IsHalo(node_index[(z*2+y)*2+x]) )
        {
          // Indicate that there is a force, probably only necessary for the unusual case of compliing without EXTERNAL_FORCES
          lbfields[node_index[(z*2+y)*2+x]].has_force = 1;

          // Add force into the lbfields structure
          double *local_f = lbfields[node_index[(z*2+y)*2+x]].force;
          
          local_f[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[0];
          local_f[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[1];
          local_f[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*delta_j[2];
        }
      }
    }
  }
}

/******************
   GetIBMInterpolatedVelocity
Very similar to the velocity interpolation done in standard Espresso, except that we add the f/2 contribution - only for CPU
*******************/

void GetIBMInterpolatedVelocity(double *p, double *const v, double *const forceAdded)
{
  index_t node_index[8], index;
  double delta[6];
  double local_rho, local_j[3], interpolated_u[3];
  double modes[19];
  int x,y,z;
  double *f;
  
  double lbboundary_mindist, distvec[3];
  double pos[3];
  
#ifdef LB_BOUNDARIES
  int boundary_no;
  int boundary_flag=-1; // 0 if more than agrid/2 away from the boundary, 1 if 0<dist<agrid/2, 2 if dist <0
  
  lbboundary_mindist_position(p, &lbboundary_mindist, distvec, &boundary_no);
  if (lbboundary_mindist>lbpar.agrid/2) {
    boundary_flag=0;
    pos[0]=p[0];
    pos[1]=p[1];
    pos[2]=p[2];
    
  } else if (lbboundary_mindist > 0 ) {
    boundary_flag=1;
    pos[0]=p[0] - distvec[0]+ distvec[0]/lbboundary_mindist*lbpar.agrid/2.;
    pos[1]=p[1] - distvec[1]+ distvec[1]/lbboundary_mindist*lbpar.agrid/2.;
    pos[2]=p[2] - distvec[2]+ distvec[2]/lbboundary_mindist*lbpar.agrid/2.;
    
  } else {
    boundary_flag=2;
    v[0]= lb_boundaries[boundary_no].velocity[0]*lbpar.agrid/lbpar.tau;
    v[1]= lb_boundaries[boundary_no].velocity[1]*lbpar.agrid/lbpar.tau;
    v[2]= lb_boundaries[boundary_no].velocity[2]*lbpar.agrid/lbpar.tau;
    return; // we can return without interpolating
  }
#else
  pos[0]=p[0];
  pos[1]=p[1];
  pos[2]=p[2];
#endif
  
  /* determine elementary lattice cell surrounding the particle
   and the relative position of the particle in this cell */
  lblattice.map_position_to_lattice(pos,node_index,delta);
  
  /* calculate fluid velocity at particle's position
   this is done by linear interpolation
   (Eq. (11) Ahlrichs and Duenweg, JCP 111(17):8225 (1999)) */
  interpolated_u[0] = interpolated_u[1] = interpolated_u[2] = 0.0 ;
  // This for the f/2 contribution to the velocity
  forceAdded[0] = forceAdded[1] = forceAdded[2] = 0;
  
  for (z=0;z<2;z++) {
    for (y=0;y<2;y++) {
      for (x=0;x<2;x++) {
        
        index = node_index[(z*2+y)*2+x];
        f = lbfields[index].force_buf;
        
        // This can be done easier withouth copying the code twice
        // We probably can even set the boundary velocity directly
#ifdef LB_BOUNDARIES
        if (lbfields[index].boundary) {
          local_rho=lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid;
          local_j[0] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[0];
          local_j[1] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[1];
          local_j[2] = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid*lb_boundaries[lbfields[index].boundary-1].velocity[2];
        } else
#endif
        {
          lb_calc_modes(index, modes);
          local_rho = lbpar.rho[0]*lbpar.agrid*lbpar.agrid*lbpar.agrid + modes[0];
          
          // Add the +f/2 contribution!!
          local_j[0] = modes[1] + f[0]/2;
          local_j[1] = modes[2] + f[1]/2;
          local_j[2] = modes[3] + f[2]/2;
          
          // Keep track of the forces that we added to the fluid
          // This is necessary for communication because this part is executed for real and ghost particles
          // Later on we sum the real and ghost contributions
          const double fExt[3] = { lbpar.ext_force[0]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau, lbpar.ext_force[1]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau, lbpar.ext_force[2]*pow(lbpar.agrid,2)*lbpar.tau*lbpar.tau };
          
          forceAdded[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*(f[0]-fExt[0])/2/(local_rho);
          forceAdded[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*(f[1]-fExt[1])/2/(local_rho);
          forceAdded[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*(f[2]-fExt[2])/2/(local_rho);
          
        }
        
        // Interpolate velocity
        interpolated_u[0] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[0]/(local_rho);
        interpolated_u[1] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[1]/(local_rho);
        interpolated_u[2] += delta[3*x+0]*delta[3*y+1]*delta[3*z+2]*local_j[2]/(local_rho) ;
        
      }
    }
  }
#ifdef LB_BOUNDARIES
  if (boundary_flag==1) {
    v[0]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[0]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[0];
    v[1]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[1]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[1];
    v[2]=lbboundary_mindist/(lbpar.agrid/2.)*interpolated_u[2]+(1-lbboundary_mindist/(lbpar.agrid/2.))*lb_boundaries[boundary_no].velocity[2];
    
  } else {
    v[0] = interpolated_u[0];
    v[1] = interpolated_u[1];
    v[2] = interpolated_u[2];
  }
#else
  v[0] = interpolated_u[0];
  v[1] = interpolated_u[1];
  v[2] = interpolated_u[2];
#endif
  v[0] *= lbpar.agrid/lbpar.tau;
  v[1] *= lbpar.agrid/lbpar.tau;
  v[2] *= lbpar.agrid/lbpar.tau;
}

/************
   IsHalo
Builds a cache structure which contains a flag for each LB node whether that node is a halo node or not
Checks for halo - only for CPU
*************/

bool IsHalo(const int indexCheck)
{
  // First call --> build cache
  if ( isHaloCache == NULL )
  {
    isHaloCache = new bool[lblattice.halo_grid_volume];
    // Assume everything is a halo and correct in the next step
    for ( int i = 0; i < lblattice.halo_grid_volume; i++)
      isHaloCache[i] = true;
    // Loop through and check where indexCheck occurs
    int index = lblattice.halo_offset;
    for (int z=1; z<=lblattice.grid[2]; z++)
    {
      for (int y=1; y<=lblattice.grid[1]; y++)
      {
        for (int x=1; x<=lblattice.grid[0]; x++)
        {
          isHaloCache[index] = false;
          ++index;
        }
        index += 2; /* skip halo region */
      }
      index += 2*lblattice.halo_grid[0]; /* skip halo region */
    }
  }
  
  // Return
  return isHaloCache[indexCheck];
}

/****************
   ParticleVelocitiesFromLB_CPU
Get particle velocities from LB and set the velocity field in the particles data structure
*****************/

void ParticleVelocitiesFromLB_CPU()
{
  // Exchange halo. This is necessary because we have done LB collide-stream
  if ( lbpar.resend_halo )
  {
    halo_communication(&update_halo_comm, (char*)**lbfluid);
    lbpar.resend_halo = 0;
  }
  
  // Loop over particles in local cells
  // Here all contributions are included: velocity, external force and particle force
  for (int c = 0; c < local_cells.n; c++)
  {
    const Cell *const cell = local_cells.cell[c];
    Particle *const p  = cell->part;
    for(int j = 0; j < cell->n; j++)
      if (p[j].p.isVirtual)
      {
        double dummy[3];
        // Get interpolated velocity and store in the force (!) field
        // for later communication (see below)
        GetIBMInterpolatedVelocity(p[j].r.p, p[j].f.f, dummy);
      }
  }
  
  // Loop over particles in ghost cells
  // Here we only add the particle forces stemming from the ghosts
  for (int c = 0; c < ghost_cells.n; c++)
  {
    const Cell *const cell = ghost_cells.cell[c];
    Particle *const p  = cell->part;
    for(int j = 0; j < cell->n; j++)
      // This criterion include the halo on the left, but excludes the halo on the right
      // Try if we have to use *1.5 on the right
      if (p[j].r.p[0] >= my_left[0]-0.5*lblattice.agrid[0]
          && p[j].r.p[0] < my_right[0]+0.5*lblattice.agrid[0]
          && p[j].r.p[1] >= my_left[1]-0.5*lblattice.agrid[1]
          && p[j].r.p[1] < my_right[1]+0.5*lblattice.agrid[1]
          && p[j].r.p[2] >= my_left[2]-0.5*lblattice.agrid[2]
          && p[j].r.p[2] < my_right[2]+0.5*lblattice.agrid[2])
      {
        if (p[j].p.isVirtual)
        {
          double dummy[3];
          double force[3]; // The force stemming from the ghost particle
          GetIBMInterpolatedVelocity(p[j].r.p, dummy, force);
          
          // Rescale and store in the force field of the particle (for communication, see below)
          p[j].f.f[0] = force[0] * lbpar.agrid/lbpar.tau;
          p[j].f.f[1] = force[1] * lbpar.agrid/lbpar.tau;
          p[j].f.f[2] = force[2] * lbpar.agrid/lbpar.tau;
        }
        else { p[j].f.f[0] = p[j].f.f[1] = p[j].f.f[2] = 0; }   // Reset, necessary because we add all forces below. Also needs to be done for the real particles!
        
      }
      else { p[j].f.f[0] = p[j].f.f[1] = p[j].f.f[2] = 0; }   // Reset, necessary because we add all forces below
  }
  
  
  
  // Now the local particles contain a velocity (stored in the force field) and the ghosts contain the rest of the velocity in their respective force fields
  // We need to add these. Since we have stored them in the force, not the velocity fields, we can use the standard force communicator and then transfer to the velocity afterwards
  // Note that this overwrites the actual force which would be a problem for real particles
  // This could be solved by keeping a backup of the local forces before this operation is attempted
  ghost_communicator(&cell_structure.collect_ghost_force_comm);
  
  // Transfer to velocity field
  for (int c = 0; c < local_cells.n; c++)
  {
    const Cell *const cell = local_cells.cell[c];
    Particle *const p  = cell->part;
    for(int j = 0; j < cell->n; j++)
      if (p[j].p.isVirtual)
      {
        p[j].m.v[0] = p[j].f.f[0];
        p[j].m.v[1] = p[j].f.f[1];
        p[j].m.v[2] = p[j].f.f[2];
      }
  }
}




#endif
