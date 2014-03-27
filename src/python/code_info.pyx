include "myconfig.pxi"

def electrostatics_defined():
    IF ELECTROSTATICS == 1:
        return True
    ELSE:
        return False
    
def dipoles_defined():
    IF DIPOLES == 1:
        return True
    ELSE:
        return False

def cuda_defined():
    IF LB_GPU == 1:
        return True
    ELSE:
        return False

