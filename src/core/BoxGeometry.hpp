#ifndef ESPRESSO_BOXGEOMETRY_HPP
#define ESPRESSO_BOXGEOMETRY_HPP

#include <config.hpp>

#include <utils/Vector.hpp>

#include <bitset>
#include <cassert>
#include <stdexcept>

struct BoxGeometry {
#ifdef PARTIAL_PERIODIC
  /** Flags for all three dimensions whether pbc are applied (default). */
  std::bitset<3> m_periodic = 0b111;

  void set_periodic(unsigned coord, bool val) { m_periodic[coord] = val; }

  constexpr bool periodic(unsigned coord) const {
    assert(coord <= 2);
    return m_periodic[coord];
  }
#else
  void set_periodic(unsigned, bool val) const {
    if (!val)
      throw std::runtime_error(
          "Disabeling periodicity needs feature 'PARTIAL_PERIODIC'.");
  }
  constexpr bool periodic(int) { return true; }
#endif

  Utils::Vector3d m_length = {1, 1, 1};

  Utils::Vector3d const &length() const { return m_length; }

  void set_length(Utils::Vector3d const &box_l) { m_length = box_l; }
};

#endif // ESPRESSO_BOXGEOMETRY_HPP
