#ifndef EXTERNAL_FIELD_DETAIL_BASE_HPP
#define EXTERNAL_FIELD_DETAIL_BASE_HPP

#include <utility>

namespace FieldCoupling {
namespace detail {
template <typename Coupling, typename Field> class Base {
protected:
  Coupling m_coupling;
  Field m_field;

public:
  template <typename CouplingRef, typename FieldRef>
  Base(CouplingRef &&coupling, FieldRef &&field)
      : m_coupling(std::forward<CouplingRef>(coupling)),
        m_field(std::forward<FieldRef>(field)) {}

  Coupling const &coupling() const { return m_coupling; }
  Field const &field() const { return m_field; }
};
}
}
#endif
