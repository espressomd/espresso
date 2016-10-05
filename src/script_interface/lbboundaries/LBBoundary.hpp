#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP

#include "ScriptInterface.hpp"
#include "core/lbboundaries/LBBoundary.hpp"
#include "core/utils/Factory.hpp"
#include "shapes/Shape.hpp"

#include <memory>

namespace ScriptInterface {
namespace LBBoundaries {
class LBBoundary : public ScriptInterfaceBase {
public:
  LBBoundary() : m_lbboundary(new ::LBBoundaries::LBBoundary()) {}

  const std::string name() const override { return  "LBBoundaries:LBBoundary"; }

  VariantMap get_parameters() const override {
    return {{"velocity", m_lbboundary->velocity()},
	    {"force", m_lbboundary->force()},
#ifdef EK_BOUNDARIES
	    {"charge_density", m_lbboundary->charge_density()},
	    {"net_charge", m_lbboundary->net_charge()},
#endif
            {"shape",
		(m_shape != nullptr) ? m_shape->id() : ScriptInterface::NOT_SET}};
  }

  ParameterMap valid_parameters() const override {
    return {{"velocity", {ParameterType::VECTOR3D, true}},
            {"force", {ParameterType::VECTOR3D, true}},
#ifdef EK_BOUNDARIES
            {"charge_density", {ParameterType::DOUBLE, true}},
	    {"net_charge", {ParameterType::DOUBLE, true}},
#endif
	    {"shape", {ParameterType::OBJECT, true}}};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    if (name == "shape") {
      std::shared_ptr<ScriptInterfaceBase> so_ptr =
          ScriptInterface::get_instance(value);

      auto shape_ptr =
          std::dynamic_pointer_cast<ScriptInterface::Shapes::Shape>(so_ptr);

      /* We are expecting a ScriptInterface::Shapes::Shape here,
         throw if not. That means the assigned object had the wrong type. */
      if (shape_ptr != nullptr) {
        m_lbboundary->set_shape(shape_ptr->shape());
        /* Store a reference */
        m_shape = shape_ptr;
      } else {
        throw std::runtime_error("shape parameter expects a Shapes::Shape");
      }
    }

    SET_PARAMETER_HELPER_VECTOR3D("velocity", m_lbboundary->velocity());
    SET_PARAMETER_HELPER_VECTOR3D("force", m_lbboundary->force());
#ifdef EK_BOUNDARIES
    SET_PARAMETER_HELPER("charge_density", m_lbboundary->charge_density());
    SET_PARAMETER_HELPER("net_charge", m_lbboundary->net_charge());
#endif
  }
    std::shared_ptr<::LBBoundaries::LBBoundary> lbboundary() {
    return m_lbboundary;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::LBBoundaries::LBBoundary> m_lbboundary;

  /* Keep a reference to the shape */
  std::shared_ptr<ScriptInterfaceBase> m_shape;

  }; // class LBBoundary
} /* namespace LBBoundary */
} /* namespace ScriptInterface */
#endif
