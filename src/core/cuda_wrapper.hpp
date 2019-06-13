#ifndef CUDA_WRAPPER_HPP
#define CUDA_WRAPPER_HPP

#define LaunchKernel(kernel, blocks, threads, mem, stream, ...)          \
  do {                                                                         \
    kernel<<<blocks, threads, mem, stream>>>(__VA_ARGS__);                     \
  } while (0)
#define DYNAMIC_SHARED(type, var) extern __shared__ type var[];
#define SYMBOL(X) X

#endif // CUDA_WRAPPER_HPP
