
#ifndef IBM_VOLUME_CONSERVATION_H
#define IBM_VOLUME_CONSERVATION_H

#ifdef IMMERSED_BOUNDARY

// This function is used to set the parameters
// Also calculates and stores the reference state
int IBM_VolumeConservation_SetParams(const int bond_type, const int softID, const double kappaV);
// For reading checkpoints.
// Idea: parameters are set in the run-continue script
//       here only check consistency of bond type
int IBM_VolumeConservation_ResetParams(const int bond_type, const double volRef);

// Init function. Calculate volumes for the first time and store as references
// Called from integrate_vv
void IBM_InitVolumeConservation();

// Calculate and apply the actual forces
// Called from forces_inline.hpp
// NB: this function is for all soft particles, volume forces cannot be calculated on a per-bond basis as usual forces
void IBM_VolumeConservation();

#endif

#endif
