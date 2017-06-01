#ifndef CORE_CELL_HPP
#define CORE_CELL_HPP

#include "particle_data.hpp"

/** A cell is a \ref ParticleList representing a particle group with
    respect to the integration algorithm.
*/
// typedef ParticleList Cell;

class Cell : public ParticleList {
public:
  void resize(size_t size) {
    realloc_particlelist(static_cast<ParticleList *>(this), this->n = size);
  }

#ifdef LEES_EDWARDS
  int myIndex[3];
#endif
};

#endif
