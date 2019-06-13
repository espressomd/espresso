#ifndef CUDA_WRAPPER_HPP
#define CUDA_WRAPPER_HPP

#define DEVICE_LAUNCH(kernel, blocks, threads, mem, stream, ...)               \
  do {                                                                         \
    kernel<<<blocks, threads, mem, stream>>>(__VA_ARGS__);                     \
  } while (0)
#define DEVICE_DYNAMIC_SHARED(type, var) extern __shared__ type var[];
#define DEVICE_SYMBOL(X) X

#endif // CUDA_WRAPPER_HPP
