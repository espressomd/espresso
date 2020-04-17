#ifndef ESPRESSO_EXCLUSIONS_HPP
#define ESPRESSO_EXCLUSIONS_HPP

#include "Particle.hpp"

#ifdef EXCLUSIONS
/** Determine if the non-bonded interactions between @p p1 and @p p2 should be
 *  calculated.
 */
inline bool do_nonbonded(Particle const &p1, Particle const &p2) {
  /* check for particle 2 in particle 1's exclusion list. The exclusion list is
   * symmetric, so this is sufficient. */
  return std::none_of(p1.exclusions().begin(), p1.exclusions().end(),
                      [&p2](int id) { return p2.p.identity == id; });
}

/** Remove exclusion from particle if possible */
void delete_exclusion(Particle *part, int part2);

/** Insert an exclusion if not already set */
void add_exclusion(Particle *part, int part2);
#endif
#endif // ESPRESSO_EXCLUSIONS_HPP
