#ifndef UTILS_INDEX_HPP
#define UTILS_INDEX_HPP

#include <array>
#include <numeric>
#include <vector>

namespace Utils {

/**
 * \brief Returns the flat index for given multidimensional indices and
 * dimensions. \param unravelled_indices a container with the multidimensional
 * indices. \param dimensions a container with the corresponding dimensions.
 * \retval the flat index
 */
template <typename T, typename U>
inline size_t ravel_index(const T &unravelled_indices, const U &dimensions) {
  const auto n_dims = unravelled_indices.size();
  if (n_dims != dimensions.size()) {
    throw std::invalid_argument(
        "Index vector and dimenions vector must have same dimensions.");
  }
  std::size_t res = unravelled_indices.back();
  std::size_t temp_prod = 1;
  for (int i = unravelled_indices.size() - 2; i >= 0; --i) {
    temp_prod *= dimensions[i + 1];
    res += unravelled_indices[i] * temp_prod;
  }
  return res;
}

/**
 * \brief Returns the unraveled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * \param dimensions an int array of length ndims containing the lengths of the
 * dimensions. (Input)
 * \param ravelled_index the flat index. (Input)
 * \retval an array with length ndims where the unflat indices are written to.
 */
template <typename T, typename U>
inline const std::vector<U> unravel_index(const T &dimensions,
                                          const U ravelled_index) {
  auto const n_dims = dimensions.size();
  std::vector<U> unravelled_indices(n_dims);
  std::size_t mul = 1;
  for (int j = n_dims - 1; j >= 0; j--) {
    unravelled_indices[j] = (ravelled_index / mul) % dimensions[j];
    mul *= dimensions[j];
  }
  return unravelled_indices;
}

} // namespace Utils

#endif
