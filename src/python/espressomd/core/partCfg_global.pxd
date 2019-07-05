from core.PartCfg cimport PartCfg
 
cdef extern from "partCfg_global.hpp":
    PartCfg & partCfg()
