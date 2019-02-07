#ifndef UTILS_INDEX_HPP
#define UTILS_INDEX_HPP

namespace Utils {

template <size_t Dims>
inline size_t ravel_index(std::vector<size_t> unravelled_indices,
                          std::array<size_t, Dims> n_bins) {
  // index calculation: using the following formula for N dimensions:
  //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k)
  size_t res = unravelled_indices.back();
  for (size_t j = 0; j < unravelled_indices.size() - 1; ++j) {
    res += unravelled_indices[j] * std::accumulate(n_bins.begin() + j + 1,
                                                   n_bins.end(), 1,
                                                   std::multiplies<size_t>());
  }
  return res;
}

/**
 * \brief Returns the unraveled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * \param len_dims an int array of length ndims containing the lengths of the
 * dimensions. (Input)
 * \param ndims int denoting the number of dimensions. (Input)
 * \param flattened_index an int denoting the flat index. (Input)
 * \param unravelled_index_out an int array with length ndims where the unflat
 * indices are written to. (Output)
 */
inline void unravel_index(const int *const len_dims, const int ndims,
                          const int flattened_index,
                          int *unravelled_index_out) {
  // idea taken from
  // http://codinghighway.com/2014/02/22/c-multi-dimensional-arrays-part-2-flattened-to-unflattened-index/
  std::vector<int> mul(ndims);
  mul[ndims - 1] = 1;
  for (int j = ndims - 2; j >= 0; j--)
    mul[j] = mul[j + 1] * len_dims[j + 1];
  for (int j = 0; j < ndims; j++)
    unravelled_index_out[j] = (flattened_index / mul[j]) % len_dims[j];
}

}

#endif
