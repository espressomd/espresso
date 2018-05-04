#include <cmath>

#include "PartCfg.hpp"

namespace Utils {

inline void check_charge_neutrality(PartCfg &partCfg) {
  double total_charge = 0.0;
  double min_charge = 0.0;
  for (auto const &particle : partCfg) {
    total_charge += particle.p.q;
    double abs_charge = std::abs(particle.p.q);
    if (min_charge == 0.0 and abs_charge != 0.0)
      min_charge = abs_charge;
    if (abs_charge < min_charge and abs_charge != 0.0)
      min_charge = particle.p.q;
  }
  if (std::abs(total_charge) / min_charge > 1e-10)
    throw std::runtime_error(
                        "The system is not charge neutral. Please "
                        "neutralize the system before adding a new actor via adding "
                        "the corresponding counterions to the system. Alternatively "
                        "you can turn off the electroneutrality check via supplying "
                        "check_neutrality=False when creating the actor. In this "
                        "case you may be simulating a non-neutral system which will "
                        "affect physical observables like e.g. the pressure, the "
                        "chemical potentials of charged species or potential "
                        "energies of the system. Since simulations of non charge "
                        "neutral systems are special please make sure you know what "
                        "you are doing.");
}

}
