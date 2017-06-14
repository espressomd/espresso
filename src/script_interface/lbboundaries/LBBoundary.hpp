#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARY_HPP

#include "ScriptInterface.hpp"
#include "auto_parameters/AutoParameters.hpp"
#include "core/lbboundaries/LBBoundary.hpp"
#include "shapes/Shape.hpp"

namespace ScriptInterface {
namespace LBBoundaries {
class LBBoundary : public AutoParameters {
public:
  LBBoundary() : m_lbboundary(new ::LBBoundaries::LBBoundary()) {
    add_parameters({{"velocity", m_lbboundary->velocity()},
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

  Variant call_method(const std::string &name, const VariantMap &) {
    if (name == "get_force") {
      return m_lbboundary->get_force();
    }
  }

  const std::string name() const override { return "LBBoundaries:LBBoundary"; }

  std::shared_ptr<::LBBoundaries::LBBoundary> lbboundary() {
    return m_lbboundary;
  }

private:
  /* The actual constraint */
  std::shared_ptr<::LBBoundaries::LBBoundary> m_lbboundary;

  /* Keep a reference to the shape */
  std::shared_ptr<Shapes::Shape> m_shape;
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
