#ifndef __SHAPE_CONTAINER_HPP
#define __SHAPE_CONTAINER_HPP

#include "Shape.hpp"
#include "ObjectManager.hpp"

namespace Shapes {
  void fill_List();
  extern  ObjectManager<Shape> List;
}
#endif
