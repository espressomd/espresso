#ifndef ESPRESSO_UTILS_FLATTEN_HPP
#define ESPRESSO_UTILS_FLATTEN_HPP

#include <iterator>
#include <type_traits>

namespace Utils {
namespace detail {
template <class Container, class OutputIterator, class = void>
struct flatten_impl {
  static OutputIterator apply(Container const &c, OutputIterator out) {
    for (auto const &e : c) {
      out = flatten_impl<decltype(e), OutputIterator>::apply(e, out);
    }

    return out;
  }
};

template <class T, class OutputIterator>
struct flatten_impl<T, OutputIterator,
                    std::enable_if_t<std::is_assignable<
                        decltype(*std::declval<OutputIterator>()), T>::value>> {
  static OutputIterator apply(T const &v, OutputIterator out) {
    *out = v;
    return ++out;
  }
};
} // namespace detail

/**
 * @brief Flatten a range of ranges.
 *
 * Copy a range of ranges to an output range by subsequently
 * copying the nested ranges to the output. Arbitrary deep
 * nesting is supported, the elements are copied into the output
 * in a depth-first fashion.
 *
 * @tparam Range A Forward Range
 * @tparam OutputIterator An OutputIterator
 * @param v Input Range
 * @param out Output iterator
 */
template <class Range, class OutputIterator>
void flatten(Range const &v, OutputIterator out) {
  detail::flatten_impl<Range, OutputIterator>::apply(v, out);
}
} // namespace Utils

#endif // ESPRESSO_FLATTEN_HPP
