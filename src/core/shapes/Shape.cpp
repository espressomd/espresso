#include "Shape.hpp"

#include "Wall.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Rhomboid.hpp"

namespace Shapes {

void initialize_factory() {
  ShapeFactory::register_new("wall", ShapeFactory::builder<Wall>);
  ShapeFactory::register_new("cylinder", ShapeFactory::builder<Cylinder>);
  ShapeFactory::register_new("sphere", ShapeFactory::builder<Sphere>);
  ShapeFactory::register_new("rhomboid", ShapeFactory::builder<Rhomboid>);
}

}
