#ifndef UTILS_RANGE_HPP
#define UTILS_RANGE_HPP

#include <utility>
#include <iterator>

namespace Utils {

/**
 * @brief Helper class that wraps two iterators into
          a range object.
*/
template <typename Iterator> class Range {
  Iterator m_begin, m_end;

public:
  using iterator = Iterator;
  using value_type = typename std::iterator_traits<Iterator>::value_type;
  using difference_type = typename std::iterator_traits<Iterator>::difference_type;

  Range(Iterator begin, Iterator end)
      : m_begin(std::forward<Iterator>(begin)),
        m_end(std::forward<Iterator>(end)) {}

  Iterator begin() { return m_begin; }
  Iterator end() { return m_end; }

  bool empty() const { return m_begin == m_end; }
  difference_type size() const {
    using std::distance;
    return distance(m_begin, m_end);
  }

  bool operator==(Range const &rhs) const {
    return (m_begin == rhs.m_begin) && (m_end == rhs.m_end);
  }
};

/**
 * @brief Return a range for a pair of iterators.
 *
 * This is a convinence function so we can template
 * argument deduction to figure out the Range type.
 */
template <typename Iterator>
Range<Iterator> make_range(Iterator begin, Iterator end) {
  return Range<Iterator>(begin, end);
}

} /* namespace Utils */

#endif
