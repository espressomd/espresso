from core.partCfg_global cimport partCfg
from libcpp cimport bool

cdef extern from "utils/checks/charge_neutrality.hpp" namespace "Utils":
    bool check_charge_neutrality[ParticleRange](ParticleRange & partCfg)
