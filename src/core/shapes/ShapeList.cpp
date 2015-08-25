#include "ShapeList.hpp"
#include "ObjectManager.hpp"
#include "Wall.hpp"

namespace Shapes {
  ObjectManager<Shape> List;

  void fill_List() {
    Factory<Shape>::Instance().register_new("wall", Factory<Shape>::builder<Wall>);
    
  }
};
