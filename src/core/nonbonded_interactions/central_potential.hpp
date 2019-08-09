#ifndef CENTRAL_POTENTIAL_HPP
#define CENTRAL_POTENTIAL_HPP

#include <tuple>
#include <unordered_map>

#include <utils/Vector.hpp>
#include <utils/hash/Cantor.hpp>
#include <utils/interaction/LennardJones.hpp>

namespace NonBondedInteractions {

template <typename Potential> struct CentralPotential : Potential {
  double shift = 0.0;
  double offset = 0.0;
  double cutoff = 0.0;
  double min = 0.0;

  double energy(double r) {
    if ((r - offset <= cutoff) and (r - offset > min)) {
      return Potential::U(r - offset) + shift;
    }
    return 0.0;
  }

  Utils::Vector3d force(double r, Utils::Vector3d const &r_ij) {
    if ((r - offset <= cutoff) and (r - offset > min)) {
      return Potential::F(r - offset, r_ij);
    }
    return {};
  }
};

using NonBondedPotentials =
    std::tuple<CentralPotential<Utils::Interaction::LennardJones>>;
using PairInteractions =
    std::unordered_map<std::pair<int, int>, NonBondedPotentials,
                       Utils::Hash::CantorPairing, Utils::Hash::CantorCompare>;
} // namespace NonBondedInteractions

#endif
