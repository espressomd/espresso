.. _GPU Acceleration with CUDA

GPU Acceleration with CUDA
**************************

.. note::
    `Feature CUDA required`


|es| is capable of GPU acceleration to speed up simulations.
Not every simulation method is parallelizable or profits from 
GPU acceleration. Refer to :ref:`Available simulation methods`
to check whether your desired method can be used on the GPU.
In order to use GPU acceleration
your GPU should at least have compute capability 2.0.

List available CUDA devices
===========================

If you want to list available CUDA devices
you should call :attr:`espressomd.cuda_init.CudaInitHandle.device_list`::

    >>> cuda_handle = espressomd.cuda_init.CudaInitHandle()
    >>> cuda_handle.device_list

Selection of CUDA device
========================

When you start ``pypresso`` your first GPU should
be selected. If you wanted
If you wanted to use the second GPU, this can be done as follows::

    >>> CUDA_handle = espressomd.CUDA_init.CUDAInitHandle()
    >>> CUDA_handle.device = 1



