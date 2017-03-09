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

  const std::string name() const override { return "LBBoundaries:LBBoundary"; }

  VariantMap get_parameters() const override {
    return {{"velocity", m_lbboundary->velocity()},
            {"force", m_lbboundary->get_force()},
#ifdef EK_BOUNDARIES
            {"charge_density", m_lbboundary->charge_density()},
            {"net_charge", m_lbboundary->net_charge()},
#endif
            {"shape", (m_shape != nullptr) ? m_shape->id() : ObjectId()}};
  }

  ParameterMap valid_parameters() const override {
    return {{"velocity", {ParameterType::DOUBLE_VECTOR, 3, true}},
            {"force", {ParameterType::DOUBLE_VECTOR, 3, true}},
#ifdef EK_BOUNDARIES
            {"charge_density", {ParameterType::DOUBLE, true}},
            {"net_charge", {ParameterType::DOUBLE, true}},
#endif
            {"shape", {ParameterType::OBJECTID, true}}};
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

    SET_PARAMETER_HELPER("velocity", m_lbboundary->velocity());
    SET_PARAMETER_HELPER("force", m_lbboundary->force());
#ifdef EK_BOUNDARIES
    if (name == "charge_density") {
      m_lbboundary->set_charge_density(boost::get<double>(value));
    }
    if (name == "net_charge") {
      m_lbboundary->set_net_charge(boost::get<double>(value));
    }
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

class LBMovingBoundary : public ScriptInterfaceBase {
public:
  LBMovingBoundary() : m_lbboundary(new ::LBBoundaries::LBMovingBoundary()) {}

  const std::string name() const override { return "LBBoundaries:LBMovingBoundary"; }

  VariantMap get_parameters() const override {
    return {{"velocity"    , m_lbboundary->velocity()    },
            {"force"       , m_lbboundary->force()       },
            {"center"      , m_lbboundary->shape().pos() },
            {"radius"      , m_lbboundary->shape().rad() },
            {"torque"      , m_lbboundary->torque()      },
            {"omega"       , m_lbboundary->omega()       },
            {"quat"        , m_lbboundary->quat()        },
            {"mass"        , m_lbboundary->mass()        },
            {"anchors"     , m_lbboundary->anchors()     },
            {"rinertia"    , m_lbboundary->rinertia()    },
            {"body_torque" , m_lbboundary->body_torque() },
            {"body_force"  , m_lbboundary->body_force()  }};
  }

  ParameterMap valid_parameters() const override {
    return {{"velocity"    , {ParameterType::DOUBLE_VECTOR, 3, true } },
            {"force"       , {ParameterType::DOUBLE_VECTOR, 3, true } },
            {"center"      , {ParameterType::DOUBLE_VECTOR, 3, true } },
            {"radius"      , {ParameterType::DOUBLE          , true } },
            {"torque"      , {ParameterType::DOUBLE_VECTOR, 3, false} },
            {"omega"       , {ParameterType::DOUBLE_VECTOR, 3, false} },
            {"quat"        , {ParameterType::DOUBLE_VECTOR, 4, false} },
            {"mass"        , {ParameterType::DOUBLE          , false} },
            {"anchors"     , {ParameterType::DOUBLE_VECTOR   , false} },
            {"rinertia"    , {ParameterType::DOUBLE_VECTOR, 3, false} },
            {"body_torque" , {ParameterType::DOUBLE_VECTOR, 3, false} },
            {"body_force"  , {ParameterType::DOUBLE_VECTOR, 3, false} }};
  }

  void set_parameter(std::string const &name, Variant const &value) override {
    SET_PARAMETER_HELPER("velocity"    , m_lbboundary->velocity()    );
    SET_PARAMETER_HELPER("force"       , m_lbboundary->force()       );
    SET_PARAMETER_HELPER("center"      , m_lbboundary->shape().pos() );
    SET_PARAMETER_HELPER("radius"      , m_lbboundary->shape().rad() );
    SET_PARAMETER_HELPER("torque"      , m_lbboundary->torque()      );
    SET_PARAMETER_HELPER("omega"       , m_lbboundary->omega()       );
    SET_PARAMETER_HELPER("quat"        , m_lbboundary->quat()        );
    SET_PARAMETER_HELPER("mass"        , m_lbboundary->mass()        );
    SET_PARAMETER_HELPER("anchors"     , m_lbboundary->anchors()     );
    SET_PARAMETER_HELPER("rinertia"    , m_lbboundary->rinertia()    );
    SET_PARAMETER_HELPER("body_torque" , m_lbboundary->body_torque() );
    SET_PARAMETER_HELPER("body_force"  , m_lbboundary->body_force()  );
  }

  std::shared_ptr<::LBBoundaries::LBMovingBoundary> lbboundary() {
    return m_lbboundary;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::LBBoundaries::LBMovingBoundary> m_lbboundary;

}; // class LBMovingBoundary
} /* namespace LBBoundaries */
} /* namespace ScriptInterface */
#endif
