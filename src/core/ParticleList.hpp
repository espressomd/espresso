#ifndef ESPRESSO_CORE_PARTICLE_LIST_HPP
#define ESPRESSO_CORE_PARTICLE_LIST_HPP

#include "Particle.hpp"

#include <utils/Span.hpp>

/** List of particles. The particle array is resized using a sophisticated
 *  (we hope) algorithm to avoid unnecessary resizes.
 *  Access using \ref realloc_particlelist, ...
 */
struct ParticleList {
  ParticleList() : part{nullptr}, n{0}, max{0} {}
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;
  /** Number of particles that fit in until a resize is needed */
  int max;

  Utils::Span<Particle> particles() { return {part, static_cast<size_t>(n)}; }
};

#endif
