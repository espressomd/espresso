#ifndef SRC_INTERACTIONS_INCLUDE_INTERACTIONS_REGISTRY_HPP
#define SRC_INTERACTIONS_INCLUDE_INTERACTIONS_REGISTRY_HPP

#include "interactions/central_potential.hpp"
#include "interactions/lennard_jones.hpp"
#include "interactions/linear.hpp"

namespace Interactions {

template <class T>
using NonBondedInteraction = boost::variant<CentralPotential<Linear<T>>,
                                            CentralPotential<LennardJones<T>>>;
template <class T>
using NonBondedInteractionsOfPair = std::vector<NonBondedInteraction<T>>;

using TypePair = std::pair<std::size_t, std::size_t>;

namespace detail {
template <class T> class Cutoff : public boost::static_visitor<T> {
public:
  template <class Inter> T operator()(Inter const &inter) const {
    return inter.cutoff();
  }
};

struct CantorPairing {
  std::size_t cantor_pairing(std::size_t k1, std::size_t k2) const {
    return (k1 + k2) / 2 * (k1 + k2 + 1) + k2;
  }

  std::size_t operator()(TypePair const &p) const {
    // make it symmetric
    if (p.second > p.first) {
      return cantor_pairing(p.first, p.second);
    }
    return cantor_pairing(p.second, p.first);
  }
};
} // namespace detail

template <class T> struct Storage {
  std::unordered_map<TypePair, NonBondedInteractionsOfPair<T>,
                     detail::CantorPairing>
      data;
  void add(TypePair const &p, NonBondedInteraction<T> const &inter) {
    data[p].push_back(inter);
  }
};

template <class StorageType, class T> struct NonBondedRegistry {
  StorageType storage;

  void add(std::size_t t1, std::size_t t2,
           NonBondedInteraction<T> const &inter) {
    storage.add(std::make_pair(t1, t2), inter);
  }

  T max_cut() const {
    return std::accumulate(
        storage.data.begin(), storage.data.end(), 0.0,
        [](auto max, auto const &pair_inters) {
          return std::accumulate(
              pair_inters.second.begin(), pair_inters.second.end(), max,
              [](auto max, auto &&inter) {
                return std::max(
                    max, boost::apply_visitor(detail::Cutoff<T>{}, inter));
              });
        });
  }
};

} // namespace Interactions

#endif // SRC_INTERACTIONS_INCLUDE_INTERACTIONS_REGISTRY_HPP
