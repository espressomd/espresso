#ifndef UTILS_INDEX_HPP
#define UTILS_INDEX_HPP

#include <iterator>
#include <numeric>

#include "Vector.hpp"

namespace Utils {

/**
 * @brief Returns the flat index for given multidimensional indices and
 * dimensions.
 * @param unravelled_indices a container with the multidimensional
 * indices.
 * @param dimensions a container with the corresponding dimensions.
 * @retval the flat index
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
 * @brief Returns the unraveled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * @param[in] dimensions_begin Iterator pointing to the begin of the container
 * with the lengths of the dimensions.
 * @param[in] dimensions_end   Iterator pointing to the end of the container
 * with the lengths of the dimensions.
 * @param[out] begin_out       Outputiterator pointing to the begin of the
 * container where the result should be written to.
 * @param[in] ravelled_index   The flat index.
 */
template <typename InputIterator, typename OutputIterator, typename T>
inline void unravel_index(InputIterator dimensions_begin,
                          InputIterator dimensions_end,
                          OutputIterator begin_out, T ravelled_index) {
  auto end_out = begin_out + std::distance(dimensions_begin, dimensions_end);
  auto rbegin_in = std::make_reverse_iterator(dimensions_begin);
  auto rend_in = std::make_reverse_iterator(dimensions_end);
  auto rend_out = std::make_reverse_iterator(end_out);
  std::size_t mul = 1;
  for (; rend_in != rbegin_in; ++rend_in, ++rend_out) {
    *rend_out = (ravelled_index / mul) % (*rend_in);
    mul *= (*rend_in);
  }
}

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

/** get the linear index from the position (@p a,@p b,@p c) in a 3D grid
 *  of dimensions @p adim.
 *
 * @return           The linear index
 * @param a , b , c  Position in 3D space
 * @param adim       Dimensions of the underlying grid
 */
inline int get_linear_index(int a, int b, int c, const Vector3i &adim) {
  assert((a >= 0) && (a < adim[0]));
  assert((b >= 0) && (b < adim[1]));
  assert((c >= 0) && (c < adim[2]));

  return (a + adim[0] * (b + adim[1] * c));
}

inline int get_linear_index(const Vector3i &ind, const Vector3i &adim) {
  return get_linear_index(ind[0], ind[1], ind[2], adim);
}

/*@}*/

} // namespace Utils

#endif
