/*
 * Copyright (C) 2019-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef UTILS_CARTESIAN_PRODUCT_HPP
#define UTILS_CARTESIAN_PRODUCT_HPP

#include <iterator>
#include <utility>

namespace Utils {
namespace detail {

template <class Body, class...> struct cart_prod_impl {
  template <class... Is> void operator()(const Body &op, Is... is) const {
    op((*is)...);
  }

  void operator()(const Body &) const { ; }
};

template <class Body, class Head, class... Tail>
struct cart_prod_impl<Body, Head, Tail...> {
  template <class... Is>
  void operator()(const Head &head, const Tail... tail, Is... is) const {
    for (auto it = std::begin(head); it != std::end(head); ++it) {
      detail::cart_prod_impl<Body, Tail...>{}(tail..., is..., it);
    }
  }
};
} // namespace detail

/**
 * @brief Call op with each element of the cartesian product set of rng.
 *
 * @param op Operation to call for each element of the product set.
 * @param rng Ranges to form the product over
 */
template <typename Body, typename... ForwardRange>
void cartesian_product(const Body &op, const ForwardRange &...rng) {
  detail::cart_prod_impl<Body, ForwardRange...>{}(rng..., op);
}
} // namespace Utils

#endif
