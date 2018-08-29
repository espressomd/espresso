#ifndef CORE_UTILS_SKIP_ITERATOR_HPP
#define CORE_UTILS_SKIP_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <cassert>
#include <iterator>

namespace Utils {

/**
 * @brief Transform iterator that skips all elements
 *        in the underlying iterator if a predicate
 *        is true.
 *
 * Example:
 *    auto a = {1, 2, ,3 ,4 ,5};
 *    it = make_skip_iterator(std::begin(a), std::end(a), [](int i) { return i %
 * 2 != 0;});
 *    Then the resulting range is {2, 4}.
 */
template <typename ForwardIterator, typename Predicate,
          typename ValueType =
              typename std::iterator_traits<ForwardIterator>::value_type>
class SkipIterator : public boost::iterator_facade<
                         SkipIterator<ForwardIterator, Predicate, ValueType>,
                         ValueType, boost::forward_traversal_tag> {
  ForwardIterator m_it, m_end;
  Predicate m_p;

public:
  SkipIterator(ForwardIterator const &it, ForwardIterator const &end,
               Predicate const &pred)
      : m_it(it), m_end(end), m_p(pred) {
    /* increment until the first element that is not skipped. */
    if ((m_it != m_end) && m_p(*m_it))
      increment();
  }

private:
  friend class boost::iterator_core_access;

  void increment() {
    /* Never increment past the end */
    if (m_it == m_end)
      return;

    do {
      ++m_it;
      /* Order in the condition is important to
         avoid dereferencing the end iterator. */
    } while ((m_it != m_end) && m_p(*m_it));
  }

  bool equal(SkipIterator const &other) const { return m_it == other.m_it; }

  ValueType &dereference() const {
    assert((m_it != m_end) && (!m_p(*m_it)));
    return *m_it;
  }
};

/**
 * @brief Use template type deduction to create a SkipIterator with
 *        the correct type.
 */
template <typename ForwardIterator, typename Predicate,
          typename SkipIter = SkipIterator<ForwardIterator, Predicate>>
SkipIter make_skip_iterator(ForwardIterator const &it,
                            ForwardIterator const &end, Predicate const &pred) {
  return SkipIter(it, end, pred);
}
} // namespace Utils

#endif
