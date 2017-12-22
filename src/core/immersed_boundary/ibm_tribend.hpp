
#ifndef IBM_TRIBEND_H
#define IBM_TRIBEND_H

#ifdef IMMERSED_BOUNDARY

// DEBUG stuff
extern double maxBendingForce, maxBendingDist, maxX;

// This function is used to set the parameters
// Also calculates and stores the reference state
int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const int ind4, const tBendingMethod method, const double kb, const bool flat);
// For reading checkpoints.
// Idea: * parameters are set in the run-continue script
//       * also reference shape is recomputed there
//       * only pass kB value here to check consistency
int IBM_Tribend_ResetParams(const int bond_type, const double kb);

// This function calculates and adds the actual force
void IBM_Tribend_CalcForce(Particle *p1, const int numPartners, Particle **const partners, const Bonded_ia_parameters &iaparams);

#endif

#endif
