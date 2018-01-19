
#ifndef IBM_TRIEL_H
#define IBM_TRIEL_H

#ifdef IMMERSED_BOUNDARY

// This function is used to set the parameters
// Also calculates and stores the reference state
int IBM_Triel_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const double max, const tElasticLaw elasticLaw, const double k1, const double k2);
// For reading checkpoints.
// Idea: * parameters are set in the run-continue script
//       * also reference shape is recomputed there
//       * only pass two values here to check consistency
int IBM_Triel_ResetParams(const int bond_type, const double k1, const double l0);

// This function calculates and adds the actual force
int IBM_Triel_CalcForce(Particle *p1,Particle *p2, Particle *p3, Bonded_ia_parameters *iaparams);

#endif

#endif
