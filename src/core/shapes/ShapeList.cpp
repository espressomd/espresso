#include "ShapeList.hpp"
#include "ObjectManager.hpp"

#include "Cylinder.hpp"
#include "HollowCone.hpp"
#include "Maze.hpp"
#include "Plane.hpp"
#include "Pore.hpp"
#include "Rhomboid.hpp"
#include "Slitpore.hpp"
#include "Sphere.hpp"
#include "SpheroCylinder.hpp"
#include "Stomatocyte.hpp"
#include "Wall.hpp"

namespace Shapes {
  ObjectManager<Shape> List;

  void init_factory() {
    Factory<Shape>::Instance().register_new("wall", Factory<Shape>::builder<Wall>);
    Factory<Shape>::Instance().register_new("cylinder", Factory<Shape>::builder<Cylinder>);
    // Factory<Shape>::Instance().register_new("hollow_cone", Factory<Shape>::builder<HollowCone>);
    // Factory<Shape>::Instance().register_new("maze", Factory<Shape>::builder<Maze>);
    // Factory<Shape>::Instance().register_new("plane", Factory<Shape>::builder<Plane>);
    // Factory<Shape>::Instance().register_new("pore", Factory<Shape>::builder<Pore>);
    // Factory<Shape>::Instance().register_new("rhomboid", Factory<Shape>::builder<Rhomboid>);
    // Factory<Shape>::Instance().register_new("slitpore", Factory<Shape>::builder<Slitpore>);
    Factory<Shape>::Instance().register_new("sphere", Factory<Shape>::builder<Sphere>);
    // Factory<Shape>::Instance().register_new("spherocylinder", Factory<Shape>::builder<SpheroCylinder>);
    // Factory<Shape>::Instance().register_new("stomatocyte", Factory<Shape>::builder<Stomatocyte>);
  }
};
