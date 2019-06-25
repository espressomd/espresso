#ifndef ESPRESSO_LOCALBOX_HPP
#define ESPRESSO_LOCALBOX_HPP

#include <utils/Vector.hpp>

template <class T> class LocalBox {
public:
  /** Dimensions of the box a single node is responsible for. */
  Utils::Vector<T, 3> local_box_l = {1, 1, 1};
  /** Left (bottom, front) corner of this nodes local box. */
  Utils::Vector<T, 3> my_left_ = {0, 0, 0};

  auto const& my_left() const { return my_left_; }
  auto const& my_right() const { return my_right_; }
  auto const& length() const { return local_box_l; }

  /** Right (top, back) corner of this nodes local box. */
  Utils::Vector<T, 3> my_right_ = {1, 1, 1};

  Utils::Vector<int, 6> boundary_ = {};
  auto const&boundary() const { return boundary_; }
};

#endif // ESPRESSO_LOCALBOX_HPP
