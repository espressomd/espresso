#ifndef UTILS_HASH_CANTOR
#define UTILS_HASH_CANTOR

#include <utility>

namespace Utils {
namespace Hash {

struct CantorPairing {
  std::size_t cantor_pairing(int k1, int k2) const {
    return static_cast<std::size_t>(0.5 * (k1 + k2) * (k1 + k2 + 1) + k2);
  }

  std::size_t operator()(std::pair<int, int> const &p) const {
    // make it symmetric
    if (p.second > p.first) {
      return cantor_pairing(p.first, p.second);
    }
    return cantor_pairing(p.second, p.first);
  }
};

struct CantorCompare {
  bool operator()(std::pair<int, int> const &lhs, std::pair<int, int> const &rhs) const {
    return (lhs.first == rhs.first and lhs.second == rhs.second) or (lhs.first == rhs.second and lhs.second == rhs.first);
  }
};

} // namespace
} // namespace

#endif
