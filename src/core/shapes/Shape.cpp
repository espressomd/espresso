#include "Shape.hpp"

#include "Wall.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Rhomboid.hpp"

namespace Shapes {

void initialize_factory() {
  ShapeFactory::Instance().register_new("wall", ShapeFactory::builder<Wall>);
  ShapeFactory::Instance().register_new("cylinder", ShapeFactory::builder<Cylinder>);
  ShapeFactory::Instance().register_new("sphere", ShapeFactory::builder<Sphere>);
  ShapeFactory::Instance().register_new("rhomboid", ShapeFactory::builder<Rhomboid>);
}

}
