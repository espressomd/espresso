cimport cellsystem
cimport global_variables

cdef class Cellsystem:
    def __init__(self):
        pass
    
    def setDomainDecomposition(self, useVerletList=True):
        useVerletList=bool(useVerletList)
        
        ### should work with global_variables.dd
        if useVerletList:
            # global_variables.dd.use_vList = 1
            pass
        else :
            # global_variables.dd.use_vList = 0
            pass
        return True
    
    def setNsquare(self):
        pass
    
    def setLayered(self, nLayers=""):
        pass
