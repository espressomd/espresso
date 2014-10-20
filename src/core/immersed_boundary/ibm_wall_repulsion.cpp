
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "communication.hpp"
#include "lb-boundaries.hpp"
#include "immersed_boundary/ibm_wall_repulsion.hpp"

/*************
   IBM_WallRepulsion_CalcForce
Calculate the repulsion and add it to the particle
 **************/

void IBM_WallRepulsion_CalcForce(Particle *p1, Bonded_ia_parameters *iaparams)
{
  if ( p1->p.isVirtual )
  {
    // Distance from next wall
    double minDist;
    double distVec[3];
    int boundaryNum;
    lbboundary_mindist_position(p1->r.p, &minDist, distVec, &boundaryNum);
    if ( minDist < 1 )
    {
      const double k = iaparams->p.ibmWallRepulsionParameters.kappaWall;
      // Use a r^(-8) to approximate an exp(-r)
      const double dist2 = minDist*minDist;
      const double dist4 = dist2*dist2;
      const double dist8 = dist4*dist4;
      const double f = k / dist8;
      
      
      // Repulsive force
      p1->f.f[0] += distVec[0] * k;
      p1->f.f[1] += distVec[1] * k;
      p1->f.f[2] += distVec[2] * k;
    }
  }
}

/****************
  IBM_WallRepulsion_ResetParams
 *****************/

int IBM_WallRepulsion_ResetParams(const int bond_type)
{
  // Check if bond exists and is of correct type
  if ( bond_type >= n_bonded_ia ) return ES_ERROR;
  if ( bonded_ia_params[bond_type].type != BONDED_IA_IBM_WALL_REPULSION ) return ES_ERROR;
  
  //Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);
  
  return ES_OK;
}

/***********
   IBM_WallRepulsion_SetParams
************/

int IBM_WallRepulsion_SetParams(const int bond_type, const double kappaWall)
{
  // Create bond
  make_bond_type_exist(bond_type);
  
  // General bond parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_WALL_REPULSION;
  bonded_ia_params[bond_type].num = 1;        // This means that Espresso requires one bond partner. Here we simply ignore it, but Espresso cannot handle 0.
  
  // Specific stuff
  bonded_ia_params[bond_type].p.ibmWallRepulsionParameters.kappaWall = kappaWall;
  
  //Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);
  
  return ES_OK;
}


#endif
