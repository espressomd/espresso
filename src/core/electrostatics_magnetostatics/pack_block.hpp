#ifndef ESPRESSO_PACK_BLOCK_HPP
#define ESPRESSO_PACK_BLOCK_HPP

#include <algorithm>

/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-block with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    in 'row-major-order' or 'C-order', that means the first index is
 *    changing slowest when running through the linear array. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index
 *    li = i2 + size[2] * (i1 + (size[1]*i0))
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
template <class T>
T *pack_block(T const *in, T *out, const int *start, const int *size,
              const int *dim, int element) {
  /* offsets for indices in input grid */
  auto const m_in_offset = element * dim[2];
  auto const s_in_offset = element * (dim[2] * (dim[1] - size[1]));
  /* Jump to start position */
  in += element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      out = std::copy_n(in, element * size[2], out);
      in += m_in_offset;
    }
    in += s_in_offset;
  }

  return out;
}

/**
 * @brief Functor that returns it's second argument.
 */
struct SecondArgument {
  template <class U, class T, class... Args>
  auto operator()(U, T &&a, Args...) const {
    return std::forward<T>(a);
  }
};

/** unpack a 3d-grid input block (size[3]) into an output 3d-grid
 *  with dimension dim[3] at start position start[3].
 *
 *  see also \ref fft_pack_block.
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 */
template <class T, class BinaryOp = SecondArgument>
const T *unpack_block(T const *in, T *out, const int *start, const int *size,
                      const int *dim, int element, BinaryOp op = {}) {
  /* copy size */
  auto const copy_size = element * size[2];
  /* offsets for indices in output grid */
  auto const m_out_offset = element * dim[2];
  auto const s_out_offset = element * (dim[2] * (dim[1] - size[1]));
  /* Jump to the start of the block */
  out += element * (start[2] + dim[2] * (start[1] + dim[1] * start[0]));

  for (int s = 0; s < size[0]; s++) {
    for (int m = 0; m < size[1]; m++) {
      std::transform(out, out + copy_size, in, out, op);
      in += copy_size;
      out += m_out_offset;
    }
    out += s_out_offset;
  }

  return in;
}

#endif // ESPRESSO_PACK_BLOCK_HPP
