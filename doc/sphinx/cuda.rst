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

Selection of Cuda device
========================

When you start ``pypresso`` your first GPU should
be selected. If you wanted
In order to use GPU acceleration you first need to select
the desired Cuda device. 
If you wanted to use the first GPU, this can be done as follows::

    >>> cuda_handle = espressomd.cuda_init.CudaInitHandle()
    >>> cuda_handle.device = 0



