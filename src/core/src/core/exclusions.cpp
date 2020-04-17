#include "exclusions.hpp"

#include <utils/contains.hpp>

#ifdef EXCLUSIONS
void add_exclusion(Particle *part, int part2) {
  if (Utils::contains(part->exclusions(), part2))
    return;

  part->exclusions().push_back(part2);
}

void delete_exclusion(Particle *part, int part2) {
  auto &el = part->exclusions();

  el.erase(std::remove(el.begin(), el.end(), part2), el.end());
}
#endif