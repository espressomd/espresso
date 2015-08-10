#ifndef __SHAPE_CONTAINER_HPP
#define __SHAPE_CONTAINER_HPP

#include "Shape.hpp"

#include <map>
#include <memory>

namespace Shapes {

#ifdef HAVE_CXX11
  typedef std::shared_ptr<Shape> pointer_type;
#else
  typedef Shape * pointer_type;
#endif
  
  class ShapeList : public std::map<int, pointer_type>
  {
  public:
    ShapeList() : m_next_id(0) {}
    int add_shape(pointer_type c) {
      insert(std::pair<int, pointer_type>(m_next_id, c));
      return m_next_id++;
    }
    void remove_shape(int i) { erase(i); }
  private:
    int m_next_id;
  };  

extern ShapeList list;

}
  
#endif
