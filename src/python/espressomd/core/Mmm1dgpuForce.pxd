include "myconfig.pxi"
from libcpp cimport bool

from core.SystemInterface cimport SystemInterface

IF ELECTROSTATICS and MMM1D_GPU:
    cdef extern from "actor/Mmm1dgpuForce.hpp":
        ctypedef float mmm1dgpu_real
        cdef cppclass Mmm1dgpuForce:
            Mmm1dgpuForce(SystemInterface & s, mmm1dgpu_real coulomb_prefactor, mmm1dgpu_real maxPWerror, mmm1dgpu_real far_switch_radius, int bessel_cutoff)
            Mmm1dgpuForce(SystemInterface & s, mmm1dgpu_real coulomb_prefactor, mmm1dgpu_real maxPWerror, mmm1dgpu_real far_switch_radius)
            Mmm1dgpuForce(SystemInterface & s, mmm1dgpu_real coulomb_prefactor, mmm1dgpu_real maxPWerror)
            void setup(SystemInterface & s)
            void tune(SystemInterface & s, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff)
            void set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff, bool manual)
            void set_params(mmm1dgpu_real _boxz, mmm1dgpu_real _coulomb_prefactor, mmm1dgpu_real _maxPWerror, mmm1dgpu_real _far_switch_radius, int _bessel_cutoff)

            unsigned int numThreads
            unsigned int numBlocks(SystemInterface & s)

            mmm1dgpu_real host_boxz
            int host_npart
            bool need_tune

            int pairs
            mmm1dgpu_real * dev_forcePairs
            mmm1dgpu_real * dev_energyBlocks

            mmm1dgpu_real coulomb_prefactor, maxPWerror, far_switch_radius
            int bessel_cutoff

            float force_benchmark(SystemInterface & s)

            void check_periodicity()
            void activate()
            void deactivate()
