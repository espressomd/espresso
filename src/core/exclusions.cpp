#include "exclusions.hpp"

#ifdef EXCLUSIONS
void add_exclusion(Particle *part, int part2) {
  for (int i = 0; i < part->el.n; i++)
    if (part->el.e[i] == part2)
      return;

  part->el.push_back(part2);
}

void delete_exclusion(Particle *part, int part2) {
  IntList &el = part->el;

  if (!el.empty()) {
    el.erase(std::remove(el.begin(), el.end(), part2), el.end());
  };
}
#endif