#ifndef LBBOUNDARIES_LBBOUNDARY_HPP
#define LBBOUNDARIES_LBBOUNDARY_HPP

#include <memory>

#include "ScriptInterface.hpp"
#include "core/utils/Factory.hpp"
#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"
#include "Vector.hpp"

namespace LBBoundaries {
  class LBBoundary : public ScriptInterfaceBase {
public:
  LBBoundary()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_velocity(Vector3d{0, 0, 0}),
        m_force(Vector3d{0, 0, 0}) {
#ifdef EK_BOUNDARIES
    m_charge_density = 0.0;
    m_net_charge = 0.0;
#endif
  }

  /* Calculate distance from the lbboundary */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }

  const std:string name() const override { return  "LBBoudaries:LBBoundary" }

  VariantMap get_parameters() const override {
    return {{"velocity", m_lbboundary->m_velocity()},
	    {"force", m_lbboundary->m_force()},
            {"shape",
		(m_shape != nullptr) ? m_shape->id() : ScriptInterface::NOT_SET}},
#ifdef EK_BOUNDARIES
	    {"charge_density", m_lbboundary->m_charge_density()},
	      {"net_charge", m_lbboundary->m_net_charge()};
#endif
  }

  ParameterMap valid_parameters() const override {
    return {{"velocity", {ParameterType::VECTOR3D, true}},
            {"force", {ParameterType::VECTOR3D, true}},
	    {"shape", {ParameterType::OBJECT, true}}},
#ifdef EK_BOUNDARIES
            {"charge_density", {ParameterType::DOUBLE, true}},
	    {"net_charge", {ParameterType::DOUBLE, true}};
#endif
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

    SET_PARAMETER_HELPER("velocity", m_lbboundary->m_velocity());
    SET_PARAMETER_HELPER("force", m_lbboundary->m_force());
#ifdef EK_BOUNDARIES
    SET_PARAMETER_HELPER("charge_density", m_lbboundary->m_charge_density());
    SET_PARAMETER_HELPER("net_charge", m_lbboundary->m_net_charge());
#endif
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  void set_velocity(Vector3d velocity) { m_velocity = velocity; }

  Shapes::Shape const &shape() const { return *m_shape; }
  Vector3d const &velocity() { return m_velocity; }
  Vector3d const &force() { return m_force; }

#ifdef EK_BOUNDARIES //TODO: ugly. Better would be a class EKBoundaries, deriving from LBBoundaries, but that requires completely different initialization infrastructure.
  void set_charge_density(float charge_density) { m_charge_density = charge_density; }
  void set_net_charge(float net_charge) { m_net_charge = net_charge; }

  float charge_density() { return m_charge_density; }
  float net_charge() { return m_net_charge; }
#endif

private:
  /** Private methods */
  /* The actual boundary */
  std::shared_ptr<::LBBoundaries::LBBoundary> m_constraint;

  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape; //TODO: I dont like this being a pointer just to get around the virtual limitations
  Vector3d m_velocity;
  Vector3d m_force;

#ifdef EK_BOUNDARIES //TODO: ugly. Better would be a class EKBoundaries, deriving from LBBoundaries, but that requires completely different initialization infrastructure.
  float m_charge_density;
  float m_net_charge;
#endif
};

} /* namespace LBBoundaries */

#endif
