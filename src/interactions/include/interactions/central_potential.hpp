#ifndef SRC_INTERACTIONS_INCLUDE_INTERACTIONS_CENTRAL_POTENTIAL_HPP
#define SRC_INTERACTIONS_INCLUDE_INTERACTIONS_CENTRAL_POTENTIAL_HPP

#include <utils/Vector.hpp>

namespace Interactions {

template <class Inter> struct CentralPotential {
  using value_type = typename Inter::value_type;

  value_type m_cutoff;
  value_type m_offset;
  value_type m_shift;
  Inter m_inter;

  CentralPotential(value_type cutoff, value_type offset, value_type shift,
                   Inter inter)
      : m_cutoff(cutoff), m_offset(offset), m_shift(shift),
        m_inter(std::move(inter)){};

  bool within_inter_range(value_type r) const {
    auto const r_eff = r - m_offset;
    return r_eff <= m_cutoff and r_eff >= 0.0;
  }

  value_type energy(value_type r) const {
    return within_inter_range(r) ? m_inter.energy(r - m_offset) + m_shift
                                 : value_type{};
  }

  Utils::Vector<value_type, 3>
  force(value_type r, Utils::Vector<value_type, 3> const &r12) const {
    return within_inter_range(r) ? m_inter.force(r - m_offset, r12)
                                 : Utils::Vector<value_type, 3>{};
  }

  value_type cutoff() const { return m_cutoff; }
};

} // namespace Interactions

#endif // SRC_INTERACTIONS_INCLUDE_INTERACTIONS_CENTRAL_POTENTIAL_HPP
