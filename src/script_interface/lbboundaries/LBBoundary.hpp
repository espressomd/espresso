#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP

#include "ScriptInterface.hpp"
#include "auto_parameters/AutoParameters.hpp"
#include "core/lbboundaries/LBBoundary.hpp"
#include "shapes/Shape.hpp"

namespace ScriptInterface {
namespace LBBoundaries {
class LBBoundary : public AutoParameters<LBBoundary> {
public:
  LBBoundary() : m_lbboundary(new ::LBBoundaries::LBBoundary()) {
    add_parameters({{"velocity",
                     [this](Variant const &value) {
                       m_lbboundary->set_velocity(get_value<Vector3d>(value));
                     },
                     [this]() {
                       return m_lbboundary->velocity();
                     }
                    },
                    {"shape",
                     [this](Variant const &value) {
                       m_shape =
                           get_value<std::shared_ptr<Shapes::Shape>>(value);

                       if (m_shape) {
                         m_lbboundary->set_shape(m_shape->shape());
                       };
                     },
                     [this]() {
                       return (m_shape != nullptr) ? m_shape->id() : ObjectId();
                     }}});
#ifdef EK_BOUNDARIES
    add_parameters(
        {{"charge_density",
          [this](Variant const &value) {
            m_lbboundary->set_charge_density(boost::get<double>(value));
          },
          [this]() { return m_lbboundary->charge_density(); },
          VariantType::DOUBLE, 0},
         {"net_charge",
          [this](Variant const &value) {
            m_lbboundary->set_net_charge(boost::get<double>(value));
          },
          [this]() { return m_lbboundary->net_charge(); }, VariantType::DOUBLE,
          0}});
#endif
  }

  Variant call_method(const std::string &name, const VariantMap &) override {
    if (name == "get_force") {
      return m_lbboundary->get_force();
    }
    return none;
  }

  std::shared_ptr<::LBBoundaries::LBBoundary> lbboundary() {
    return m_lbboundary;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::LBBoundaries::LBBoundary> m_lbboundary;

  /* Keep a reference to the shape */
  std::shared_ptr<Shapes::Shape> m_shape;
}; // class LBBoundary
} /* namespace LBBoundary */
} /* namespace ScriptInterface */
#endif
