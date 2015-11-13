#include "Shape.hpp"

#include "Wall.hpp"
#include "Cylinder.hpp"

namespace Shapes {

void initialize_factory() {
  ShapeFactory::Instance().register_new("wall", ShapeFactory::builder<Shapes::Wall>);
  ShapeFactory::Instance().register_new("cylinder", ShapeFactory::builder<Shapes::Cylinder>);
}

}
