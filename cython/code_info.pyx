
cdef extern from *:
    void IFDEF_ELECTROSTATICS "#ifdef ELECTROSTATICS //" () 
    void ENDIF "#endif //" () 
    void ELS "#else //" () 

def electrostatics():
    IFDEF_ELECTROSTATICS()
    return True
    ELS()
    return False
    ENDIF()



