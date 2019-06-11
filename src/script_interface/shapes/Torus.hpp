#ifndef SCRIPT_INTERFACE_TORUS_WALL_HPP
#define SCRIPT_INTERFACE_TORUS_WALL_HPP

#include "Shape.hpp"
#include "core/shapes/Torus.hpp"

namespace ScriptInterface {
namespace Shapes {

class Torus : public Shape {
  using CoreShape = ::Shapes::Torus;
  std::shared_ptr<::Shapes::Torus> m_torus;

public:
  Torus() : m_torus(new ::Shapes::Torus()) {
    add_parameters(
        {{"radius", m_torus, &CoreShape::set_radius, &CoreShape::radius},
         {"tube_radius", m_torus, &CoreShape::set_tube_radius,
          &CoreShape::tube_radius},
         {"normal", m_torus, &CoreShape::set_normal, &CoreShape::normal},
         {"center", m_torus, &CoreShape::center},
         {"direction", m_torus, &CoreShape::direction}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override { return m_torus; }
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
