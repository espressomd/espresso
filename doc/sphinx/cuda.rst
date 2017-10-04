.. _GPU Acceleration with CUDA:

GPU Acceleration with CUDA
**************************

.. note::
    `Feature CUDA required`


|es| is capable of GPU acceleration to speed up simulations.
Not every simulation method is parallelizable or profits from 
GPU acceleration. Refer to :ref:`Available simulation methods`
to check whether your desired method can be used on the GPU.
In order to use GPU acceleration you need a NVIDIA GPU
and it needs to have at least compute capability 2.0.

For more information please check :attr:`espressomd._system.cu`
or :class:`espressomd.cuda_init.CudaInitHandle`.

.. _List available CUDA devices:

List available CUDA devices
===========================

If you want to list available CUDA devices
you should access :attr:`espressomd._system.cu.device_list`::

    espressomd._system.cu.device_list

This attribute is read only and will return a dictionary containing
the device id as key and the device name as its' value.

.. _Selection of CUDA device:

Selection of CUDA device
========================

When you start ``pypresso`` your first GPU should
be selected. 
If you wanted to use the second GPU, this can be done 
by setting :attr:`espressomd._system.cu.device` as follows::

    espressomd._system.cu.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.


