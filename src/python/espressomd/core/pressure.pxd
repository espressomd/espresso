from core.statistics cimport Observable_stat, Observable_stat_non_bonded

cdef extern from "pressure.hpp":
    cdef Observable_stat total_pressure
    cdef Observable_stat_non_bonded total_pressure_non_bonded
    cdef Observable_stat total_p_tensor
    cdef Observable_stat_non_bonded total_p_tensor_non_bonded
    cdef void update_pressure(int)
