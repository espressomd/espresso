#ifndef UTILS_CARTESIAN_PRODUCT_HPP
#define UTILS_CARTESIAN_PRODUCT_HPP

#include <utility>

namespace Utils {
namespace detail {

    template <class Body, class...> struct cart_prod_impl {
        template <class... Is>
        void operator()(const Body &op, Is... is) const {
                op((*is)...);
        }

        void operator()(const Body&) const {
            ;
        }
    };

template <class Body, class Head, class... Tail> struct cart_prod_impl<Body, Head, Tail...> {
  template <class... Is>
  void operator()(const Head &head, const Tail... tail, Is... is) const {
    using std::begin;
    using std::end;

    for (auto it = begin(head); it != end(head); ++it) {
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
void cartesian_product(const Body &op, const ForwardRange&... rng) {
  detail::cart_prod_impl<Body, ForwardRange...>{}(rng..., op);
}
} // namespace Utils

#endif
