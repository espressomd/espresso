#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <boost/range/algorithm/mismatch.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/numeric.hpp>

namespace Utils {

/**
 * @brief Heap-allocated tensor class.
 * @tparam T data type
 */
template <typename T> class Tensor {
public:
  using value_type = T;
  using reference = T &;
  using const_reference = T const &;
  using pointer = T *;
  using const_pointer = const value_type *;
  using iterator = pointer;
  using const_iterator = const_pointer;

  /**
   * @brief Creates a tensor with given extents
   * @param extents Size for each dimension
   */
  explicit Tensor(std::initializer_list<std::size_t> const &extents)
      : Tensor(std::begin(extents), std::end(extents)) {}

  template <typename Container>
  Tensor(Container const &extents)
      : Tensor(std::begin(extents), std::end(extents)) {}

  template <typename Iter>
  Tensor(Iter extents_begin, Iter extents_end)
      : m_extents(extents_begin, extents_end) {
    m_data = std::vector<T>(std::accumulate(m_extents.begin(), m_extents.end(),
                                            1, std::multiplies<std::size_t>()));

    boost::transform(m_extents, std::back_inserter(m_strides),
                     [stride = std::size_t(1)](std::size_t dim) mutable {
                       auto const old = stride;
                       stride *= dim;

                       return old;
                     });
  }

  /**
   * @brief Element access.
   * @param indices index for each dimension
   */
  template <typename... size_ts>
  const_reference operator()(size_ts... indices) const {
    return operator()({static_cast<std::size_t>(indices)...});
  }

  template <typename... size_ts> reference operator()(size_ts... indices) {
    return operator()({static_cast<std::size_t>(indices)...});
  }

  template <typename Container>
  const_reference operator()(Container const &indices) const {
    assert(valid_indices(indices));
    return m_data[calc_linear_index(indices)];
  }

  template <typename Container> reference operator()(Container const &indices) {
    assert(valid_indices(indices));
    return m_data[calc_linear_index(indices)];
  }

  const_reference
  operator()(std::initializer_list<std::size_t> const &indices) const {
    assert(valid_indices(indices));
    return m_data[calc_linear_index(indices)];
  }

  reference operator()(std::initializer_list<std::size_t> const &indices) {
    assert(valid_indices(indices));
    return m_data[calc_linear_index(indices)];
  }

  /**
   * @brief Bounds checked element access.
   * @param indices index for each dimension
   */
  template <typename... size_ts> const_reference at(size_ts... indices) const {
    return at({static_cast<std::size_t>(indices)...});
  }

  template <typename... size_ts> reference at(size_ts... indices) {
    return at({static_cast<std::size_t>(indices)...});
  }

  template <typename Container>
  const_reference at(Container const &indices) const {
    validate_indices(indices);
    return operator()(indices);
  }

  template <typename Container> reference at(Container const &indices) {
    validate_indices(indices);
    return operator()(indices);
  }

  const_reference at(std::initializer_list<std::size_t> const &indices) const {
    validate_indices(indices);
    return operator()(indices);
  }

  reference at(std::initializer_list<std::size_t> const &indices) {
    validate_indices(indices);
    return operator()(indices);
  }

  reference front() { return *begin(); }

  const_reference front() const { return *cbegin(); }

  reference back() { return m_data.back(); }

  const_reference back() const { return m_data.back(); }

  iterator begin() noexcept { return &(*m_data.begin()); };

  const_iterator begin() const noexcept { return &(*m_data.begin()); };

  const_iterator cbegin() const noexcept { return &(*m_data.begin()); };

  iterator end() noexcept { return begin() + m_data.size(); };

  const_iterator end() const noexcept { return begin() + m_data.size(); };

  const_iterator cend() const noexcept { return begin() + m_data.size(); };

  pointer data() noexcept { return m_data.data(); }

  const_pointer data() const noexcept { return m_data.data(); }

  auto rank() const noexcept { return m_extents.size(); }

  auto const &extents() const noexcept { return m_extents; }

private:
  /**
   * @brief Check that indices are valid.
   * The number of given indices as well as being within the extents is checked.
   * @param indices Indices to be checked.
   * @returns boolean true if valid, false else.
   * @tparam Type of container holding the indices.
   */
  template <typename Container>
  bool valid_indices(Container const &indices) const {
    if (indices.size() != m_extents.size()) {
      return false;
    }
    // check if any index exceeds an extend.
    auto const res =
        boost::mismatch(indices, m_extents, std::less<std::size_t>());
    return res.first == indices.end();
  }

  /**
   * @brief Throw if invalid indices.
   * @param indices Indices to be checked.
   * @tparam Type of container holding the indices.
   * @throws std::runtime_error if indices are invalid.
   */
  template <typename Container>
  void validate_indices(Container const &indices) const {
    if (not valid_indices(indices)) {
      throw std::runtime_error("Invalid indices");
    }
  }

  template <typename Container>
  std::size_t calc_linear_index(Container const &indices) {
    return boost::inner_product(m_strides, indices, 0);
  }

  std::vector<T> m_data;
  std::vector<std::size_t> m_extents;
  std::vector<std::size_t> m_strides;
};

} // namespace Utils
#endif
