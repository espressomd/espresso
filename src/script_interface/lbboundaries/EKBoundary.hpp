#ifndef SCRIPT_INTERFACE_EKBOUNDARIES_EKBOUNDARY_HPP
#define SCRIPT_INTERFACE_EKBOUNDARIES_EKBOUNDARY_HPP
#include "config.hpp"

#include "core/communication.hpp"
#include "core/grid_based_algorithms/lbboundaries/EKBoundary.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/shapes/Shape.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace EKBoundaries {
class EKBoundary : public AutoParameters<EKBoundary> {
public:
  EKBoundary() : m_ekboundary(new ::EKBoundaries::EKBoundary()) {
    add_parameters({{"shape",
                     [this](Variant const &value) {
                       m_shape =
                           get_value<std::shared_ptr<Shapes::Shape>>(value);

                       if (m_shape) {
                         m_ekboundary->set_shape(m_shape->shape());
                       };
                     },
                     [this]() { return m_shape; }}});
  }

  std::shared_ptr<::EKBoundaries::EKBoundary> ekboundary() {
    return m_ekboundary;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::EKBoundaries::EKBoundary> m_ekboundary;

  /* Keep a reference to the shape */
  std::shared_ptr<Shapes::Shape> m_shape;
};
} // namespace EKBoundaries
} /* namespace ScriptInterface */

#endif // SCRIPT_INTERFACE_EKBOUNDARIES_EKBOUNDARY_HPP
