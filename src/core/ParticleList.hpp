#ifndef ESPRESSO_CORE_PARTICLE_LIST_HPP
#define ESPRESSO_CORE_PARTICLE_LIST_HPP

#include "Particle.hpp"

#include <utils/Span.hpp>
#include <utils/memory.hpp>

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

  int resize(int size) {
/** granularity of the particle buffers in particles */
    constexpr int PART_INCREMENT  = 8;

    assert(size >= 0);
    int old_max = max;
    Particle *old_start = part;

    if (size < max) {
      if (size == 0)
        /* to be able to free an array again */
        max = 0;
      else
        /* shrink not as fast, just lose half, rounded up */
        max =
            PART_INCREMENT *
            (((max + size + 1) / 2 + PART_INCREMENT - 1) / PART_INCREMENT);
    } else
      /* round up */
      max = PART_INCREMENT * ((size + PART_INCREMENT - 1) / PART_INCREMENT);
    if (max != old_max)
      part = Utils::realloc(part, sizeof(Particle) * max);
    return part != old_start;
  }
};

/** Allocate storage for local particles and ghosts. This version
    does \em not care for the bond information to be freed if necessary.
    \param plist the list on which to operate
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT.
    \return true iff particle addresses have changed */
inline int realloc_particlelist(ParticleList *l, int size) {
  return assert(l), l->resize(size);
}

#endif
