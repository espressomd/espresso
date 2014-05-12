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
        
        # grid.h::node_grid
        mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC)
        
        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True
    
    def setNsquare(self):
        mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE)
        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True
    
    def setLayered(self, nLayers=""):
        if nLayers!="":
            if not isinstance(nLayers, int):
                raise ValueError("layer height should be positive")
            global n_layers
            n_layers=int(nLayers)
            global determine_n_layers
            determine_n_layers=0
        
        if (node_grid[0] != 1 or node_grid[1] != 1):
            node_grid[0] = node_grid[1] = 1
            node_grid[2] = n_nodes;
            err = mpi_bcast_parameter(FIELD_NODEGRID)
        else:
            err = 0

        if not err:
            mpi_bcast_cell_structure(CELL_STRUCTURE_LAYERED)
        
        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True
            
