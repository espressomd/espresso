#ifndef ESPRESSO_CUDA_HOST_ALLOCATOR_HPP
#define ESPRESSO_CUDA_HOST_ALLOCATOR_HPP

#include <type_traits>
#include <vector>

void cuda_malloc_host(void **p, std::size_t n);
void cuda_free_host(void *p);

/**
 * @brief Allocator that uses CUDA to allocate CPU memory.
 *
 * Using the CUDA allocator can have performance benefits,
 * because it returns pinned memory that is suitable for
 * DMA.
 *
 * @tparam T Type to allocate memory for.
 */
template <class T> struct CudaHostAllocator {
  using value_type = T;
  using pointer = T *;
  using reference = T &;
  using const_reference = std::add_const_t<reference>;

  CudaHostAllocator() noexcept = default;
  template <class U> explicit CudaHostAllocator(const CudaHostAllocator<U> &) {}
  template <class U> bool operator==(const CudaHostAllocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const CudaHostAllocator<U> &) const {
    return false;
  }

  T *allocate(const std::size_t n) const {
    T *result{};

    cuda_malloc_host(reinterpret_cast<void **>(&result),
                     n * sizeof(value_type));

    return result;
  }
  void deallocate(T *const p, std::size_t) const noexcept { cuda_free_host(p); }
};

template <class T> using pinned_vector = std::vector<T, CudaHostAllocator<T>>;
#endif
