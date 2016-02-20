#include "Shape.hpp"

#include "Wall.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Rhomboid.hpp"

namespace Shapes {

void initialize_factory() {
  Factory::register_new("wall", Factory::builder<Wall>);
  Factory::register_new("cylinder", Factory::builder<Cylinder>);
  Factory::register_new("sphere", Factory::builder<Sphere>);
  Factory::register_new("rhomboid", Factory::builder<Rhomboid>);
}

}
