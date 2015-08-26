#ifndef __SHAPE_CONTAINER_HPP
#define __SHAPE_CONTAINER_HPP

#include "Shape.hpp"
#include "ObjectManager.hpp"

namespace Shapes {
  void init_factory();
  extern  ObjectManager<Shape> List;
}
#endif
