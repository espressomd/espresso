#ifndef ESPRESSO_PACK_BLOCK_HPP
#define ESPRESSO_PACK_BLOCK_HPP

#include <utils/index.hpp>

#include <algorithm>

namespace detail {
template <size_t slow_dim, size_t mid_dim, size_t fast_dim, class T>
T *pack_block(T const *in, T *out, const int *start, const int *size,
              const int *dim, int element) {
  /* offsets for indices in input grid */
  auto const m_in_offset = element * dim[fast_dim];
  auto const s_in_offset =
      element * (dim[fast_dim] * (dim[mid_dim] - size[mid_dim]));
  /* Jump to start position */
  in += element *
        (start[fast_dim] +
         dim[fast_dim] * (start[mid_dim] + dim[mid_dim] * start[slow_dim]));
  /* Copy size */
  auto const copy_size = element * size[fast_dim];

  for (int s = 0; s < size[slow_dim]; s++) {
    for (int m = 0; m < size[mid_dim]; m++) {
      out = std::copy_n(in, copy_size, out);
      in += m_in_offset;
    }
    in += s_in_offset;
  }

  return out;
}

template <size_t slow_dim, size_t mid_dim, size_t fast_dim, class T,
          class BinaryOp>
const T *unpack_block(T const *in, T *out, const int *start, const int *size,
                      const int *dim, int element, BinaryOp op) {
  /* copy size */
  auto const copy_size = element * size[fast_dim];
  /* offsets for indices in output grid */
  auto const m_out_offset = element * dim[fast_dim];
  auto const s_out_offset =
      element * (dim[fast_dim] * (dim[mid_dim] - size[mid_dim]));
  /* Jump to the start of the block */
  out += element *
         (start[fast_dim] +
          dim[fast_dim] * (start[mid_dim] + dim[mid_dim] * start[slow_dim]));

  for (int s = 0; s < size[slow_dim]; s++) {
    for (int m = 0; m < size[mid_dim]; m++) {
      std::transform(out, out + copy_size, in, out, op);
      in += copy_size;
      out += m_out_offset;
    }
    out += s_out_offset;
  }

  return in;
}
} // namespace detail

/** pack a block (size[3] starting at start[3]) of an input 3d-grid
 *  with dimension dim[3] into an output 3d-block with dimension size[3].
 *
 *    The block with dimensions (size[0], size[1], size[2]) is stored
 *    according to the specified memory order. The
 *    element (i0 (slow), i1 (mid), i2 (fast)) has the linear index
 *    li = i2 + size[2] * (i1 + (size[1]*i0))
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 *  \param[in]  memory_order column- or row-major
 */
template <class T>
T *pack_block(T const *in, T *out, const int *start, const int *size,
              const int *dim, int element, Utils::MemoryOrder memory_order) {
  switch (memory_order) {
  case Utils::MemoryOrder::ROW_MAJOR:
    return detail::pack_block<0, 1, 2>(in, out, start, size, dim, element);
  case Utils::MemoryOrder::COLUMN_MAJOR:
    return detail::pack_block<2, 1, 0>(in, out, start, size, dim, element);
  default:
    throw std::runtime_error("Unknown memory order.");
  }
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
 *  see also \ref pack_block.
 *
 *  \param[in]  in      pointer to input 3d-grid.
 *  \param[out] out     pointer to output 3d-grid (block).
 *  \param[in]  start   start index of the block in the in-grid.
 *  \param[in]  size    size of the block (=dimension of the out-grid).
 *  \param[in]  dim     size of the in-grid.
 *  \param[in]  element size of a grid element (e.g. 1 for Real, 2 for Complex).
 *  \param[in]  memory_order column- or row-major
 */
template <class T, class BinaryOp = SecondArgument>
const T *unpack_block(T const *in, T *out, const int *start, const int *size,
                      const int *dim, int element,
                      Utils::MemoryOrder memory_order, BinaryOp &&op = {}) {
  switch (memory_order) {
  case Utils::MemoryOrder::ROW_MAJOR:
    return detail::unpack_block<0, 1, 2>(in, out, start, size, dim, element,
                                         std::forward<BinaryOp>(op));
  case Utils::MemoryOrder::COLUMN_MAJOR:
    return detail::unpack_block<2, 1, 0>(in, out, start, size, dim, element,
                                         std::forward<BinaryOp>(op));
  default:
    throw std::runtime_error("Unknown memory order.");
  }
}

#endif // ESPRESSO_PACK_BLOCK_HPP
