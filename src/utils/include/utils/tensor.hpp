#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <algorithm>
#include <array>
#include <numeric>
#include <vector>

#include <boost/range/algorithm/transform.hpp>
#include <boost/range/numeric.hpp>

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
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  /**
   * @brief Creates a tensor with given extends
   * @param extends Size for each dimension
   */
  explicit Tensor(std::initializer_list<std::size_t> const &extends)
      : m_extends(extends.begin(), extends.end()) {
    m_data = std::vector<T>(std::accumulate(extends.begin(), extends.end(), 1,
                                            std::multiplies<std::size_t>()));

    boost::transform(extends, std::back_inserter(m_strides),
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

  const_reference
  operator()(std::initializer_list<std::size_t> const &indices) const {
    validate_indices(indices);
    return m_data[calc_linear_index(indices)];
  }

  reference operator()(std::initializer_list<std::size_t> const &indices) {
    validate_indices(indices);
    return m_data[calc_linear_index(indices)];
  }

  reference front() { return *begin(); }

  const_reference front() const { return *cbegin(); }

  reference back() { return *(end() - 1); }

  const_reference back() const { return *(cend() - 1); }

  iterator begin() noexcept { return m_data.begin(); };

  const_iterator begin() const noexcept { return m_data.begin(); };

  const_iterator cbegin() const noexcept { return m_data.begin(); };

  iterator end() noexcept { return m_data.end(); };

  const_iterator end() const noexcept { return m_data.end(); };

  const_iterator cend() const noexcept { return m_data.end(); };

  pointer data() noexcept { return m_data.data(); }

  const_pointer data() const noexcept { return m_data.data(); }

  auto rank() const noexcept { return m_extends.size(); }

  auto extends() const noexcept { return m_extends; }

private:
  void
  validate_indices(std::initializer_list<std::size_t> const &indices) const {
    if (indices.size() != m_extends.size()) {
      throw std::runtime_error("Number of indices have to match rank");
    }
    // check if any index exceeds an extend.
    auto const res =
        std::mismatch(indices.begin(), indices.end(), m_extends.begin(),
                      [](std::size_t a, std::size_t b) { return a < b; });
    if (not(res.first == indices.end())) {
      throw std::out_of_range("invalid index");
    }
  }

  std::size_t
  calc_linear_index(std::initializer_list<std::size_t> const &indices) {
    return boost::inner_product(m_strides, indices, 0);
  }

  std::vector<T> m_data;
  std::vector<std::size_t> m_extends;
  std::vector<std::size_t> m_strides;
};

#endif
