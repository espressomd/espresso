#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <array>
#include <numeric>
#include <vector>

#include <boost/range/algorithm/transform.hpp>
#include <boost/range/numeric.hpp>

/**
 * @brief Heap-allocated tensor class.
 * @tparam T data type
 * @tparam t_rank number of dimensions of the tensor
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
   * @brief Creates a tensor of given rank and sizes
   * @param dimensions Size for each dimension
   * @param rank Rank of the tensor
   */
  explicit Tensor(std::initializer_list<std::size_t> const &dimensions)
      : m_rank(dimensions.size()),
        m_dimensions(dimensions.begin(), dimensions.end()) {
    m_data =
        std::vector<T>(std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                       std::multiplies<std::size_t>()));

    boost::transform(dimensions, std::back_inserter(m_strides),
                     [stride = 1](size_t dim) mutable {
                       auto const old = stride;
                       stride *= dim;

                       return old;
                     });
  }

  /**
   * @brief Element access.
   * @param indices index for each dimension
   */
  const_reference
  operator()(std::initializer_list<std::size_t> const &indices) const {
    auto const index = boost::inner_product(m_strides, indices, 0);
    assert(index < m_data.size());
    return m_data[index];
  }

  reference operator()(std::initializer_list<std::size_t> const &indices) {
    auto const index = boost::inner_product(m_strides, indices, 0);
    assert(index < m_data.size());
    return m_data[index];
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

  auto rank() const noexcept { return m_rank; }

  auto dimensions() const noexcept { return m_dimensions; }

private:
  std::size_t m_rank;
  std::vector<T> m_data;
  std::vector<std::size_t> m_dimensions;
  std::vector<std::size_t> m_strides;
};

#endif
