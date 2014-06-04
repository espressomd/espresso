cimport thermostat

cdef class Thermostat:
    def __init__(self):
        pass

    def getStatus(self):
        """Returns the thermostat status. Equivalent to the tcl command 'thermostat' with no arguments"""
        
        if temperature == -1:
            return "{ not initialized }"
        if thermo_switch == THERMO_OFF:
            return "{ off }"
        if thermo_switch and THERMO_LANGEVIN :
            return "{ langevin "+str(temperature)+" "+str(langevin_gamma)+" }"
    
    def turnOff(self):
        """Turns off all the thermostat and sets all the thermostat variables to zero"""
        
        global temperature
        temperature=0.
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        global langevin_gamma
        langevin_gamma=0.
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        global thermo_switch
        thermo_switch = THERMO_OFF
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        # here other thermostats stuff
        return True
    
    def setLangevin(self, _temperature="", _gamma=""):
        """Sets the Langevin thermostat with required parameters 'temperature' 'gamma'"""
        
        if _temperature=="" or _gamma=="":
            raise ValueError("wrong # args:  should be\n\"thermostat langevin <temp> <gamma>\"")
        if not isinstance(_temperature, float) or not isinstance(_gamma, float) or float(_temperature)<0. or float(_gamma)<0.:
            raise ValueError("temperature and gamma must be positive numbers")
        global temperature
        temperature=float(_temperature)
        global langevin_gamma
        langevin_gamma=float(_gamma)
        global thermo_switch
        thermo_switch = ( thermo_switch or THERMO_LANGEVIN )
        mpi_bcast_parameter(FIELD_THERMO_SWITCH)
        mpi_bcast_parameter(FIELD_TEMPERATURE)
        mpi_bcast_parameter(FIELD_LANGEVIN_GAMMA)
        return True
