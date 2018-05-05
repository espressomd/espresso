#include <cmath>
#include <limits>

namespace Utils {

template<typename ParticleRange, typename Particle>
inline void check_charge_neutrality(ParticleRange &prange) {
  std::vector<double> charges;
  std::transform(prange.begin(), prange.end(), std::back_inserter(charges),
          [](Particle const &p){return p.p.q;});
  double total_charge = std::accumulate(charges.begin(), charges.end(), 0.0);
  std::sort(charges.begin(), charges.end(),
          [](double a, double b) {return std::abs(a) < std::abs(b);});
  //Find first non-zero element.
  auto result = std::find_if(charges.begin(), charges.end(), [](double a){return std::abs(a) > std::numeric_limits<double>::epsilon();});
  if (result != charges.end()) {
    double min_abs_non_zero_charge = std::abs(*result);
    if (std::abs(total_charge) / min_abs_non_zero_charge > std::numeric_limits<double>::epsilon())
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
}
