#ifndef ESPRESSO_LOCALBOX_HPP
#define ESPRESSO_LOCALBOX_HPP

#include <utils/Vector.hpp>

template <class T> class LocalBox {
public:
  /** Dimensions of the box a single node is responsible for. */
  Utils::Vector<T, 3> local_box_l = {1, 1, 1};
  /** Left (bottom, front) corner of this nodes local box. */
  Utils::Vector<T, 3> my_left = {0, 0, 0};
  /** Right (top, back) corner of this nodes local box. */
  Utils::Vector<T, 3> my_right = {1, 1, 1};
};

#endif // ESPRESSO_LOCALBOX_HPP
