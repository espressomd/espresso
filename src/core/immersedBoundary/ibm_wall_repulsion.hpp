
#ifndef IBM_WALL_REPULSION_H
#define IBM_WALL_REPULSION_H

#ifdef IMMERSED_BOUNDARY

// These functions are used to set the parameters
// Second one is for reading checkpoints
int IBM_WallRepulsion_SetParams(const int bond_type, const double kappaWall);
int IBM_WallRepulsion_ResetParams(const int bond_type);

// This function calculates and adds the actual force
void IBM_WallRepulsion_CalcForce(Particle *p_ind1, Bonded_ia_parameters *iaparams);

#endif

#endif
