include "myconfig.pxi"
import lb
import numpy as np

IF ELECTROKINETICS:

    cdef class Electrokinetics(lb.HydrodynamicInteraction):
        
        def validate_params(self):
            default_params = self.default_params()
            
            if not (self._params["stencil"] in ["linkcentered", "nonlinear", "nodecentered"]):
                raise ValueError("stencil has to be 'linkcentered', 'nonlinear' or 'nodecentered'.")

            if not (self._params["fluid_coupling"] in ["friction", "estatics"]):
                raise ValueError("fluid_coupling has to be 'friction' or 'estatics'.")
                
        
        def valid_keys(self):
            return "agrid", "lb_density", "viscosity", "friction", "bulk_viscosity", "gamma_even", "gamma_odd", "T", "bjerrum_length", "stencil", "advection", "fluid_coupling"

        def required_keys(self):
            return ["agrid", "lb_density", "viscosity", "friction", "T", "bjerrum length"]

        def default_params(self):
            return {"agrid": -1,
                    "lb_density": -1,
                    "viscosity": -1,
                    "bulk_viscosity": -1,
                    "gamma_odd": 0.0,
                    "gamma_even": 0.0,
                    "friction": 0.0,
                    "T": -1,
                    "bjerrum_length": -1,
                    "stencil": 0,
                    "advection": True,
                    "fluid_coupling": True}

        def _get_params_from_es_core(self):
            if ek_parameters.stencil == 0:
                stencil = "linkcentered"
            elif ek_parameters.stencil == 1:
                stencil = "nonlinear"
            elif ek_parameters.stencil == 2:
                stencil = "nodecentered"
            else: 
                raise Exception("Value of stencil could not be identified.")
                
            if ek_parameters.fluidcoupling_ideal_contribution == True:
                fluid_coupling = "friction"
            else:
                fluid_coupling = "estatics"

            return {"agrid": ek_parameters.agrid,
                    "lb_density": ek_parameters.lb_density,
                    "viscosity": ek_parameters.viscosity,
                    "bulk_viscosity": ek_parameters.bulk_viscosity,
                    "gamma_odd": ek_parameters.gamma_odd,
                    "gamma_even": ek_parameters.gamma_even,
                    "friction": ek_parameters.friction,
                    "T": ek_parameters.T,
                    "bjerrum_length":ek_parameters.bjerrumlength,
                    "stencil": stencil,
                    "advection": ek_parameters.advection,
                    "fluid_coupling": fluid_coupling}

            
        def _set_params_in_es_core(self):
            if self._params["stencil"] = "linkcentered":
                ek_set_stencil(0)
            elif self._params["stencil"] = "nonlinear":
                ek_set_stencil(1)
            elif self._params["stencil"] = "nodecentered":
                ek_set_stencil(2)

            if self._params["fluid_coupling"] = "friction":
                ek_set_fluidcoupling(True)
            elif self._params["fluid_coupling"] = "estatics":
                ek_set_fluidcoupling(False)

            ek_set_agrid(self._params["agrid"])
            ek_set_lb_density(self._parms["lb_density"])
            ek_set_viscosity(self._params["viscosity"])
            ek_set_friction(self._params["friction"])
            ek_set_T(self._params["T"])
            ek_set_bjerrumlength(self._params["bjerrum_length"])
            ek_set_bulk_viscosity(self._params["bulk_viscosity"])
            ek_set_gamma_odd(self._params["gamma_odd"])
            ek_set_gamma_even(self._params["gamma_even"])
            ek_set_advection(self._params["advection"])

            
        def set_density(self, species=None, density=None, node=None):
            if spicies == None or density == None:
                raise ValueError("species and density has to be set.")
            if not isinstance(int, species):
                raise ValueError("species needs to be an integer.")
            if node == None:
                ek_set_density(species, density)
            else:
                if not (isinstance(list, node) or isinstance(np.ndarray, node)):
                    if len(coords) != 3:
                        raise ValueError("coords has to be an array of length 3 of integers.")
                ek_node_set_density(species, node[0], node[1], node[2], density)

                    
        def _activate_method(self):
            self._set_params_in_es_core()
            err = ek_init()
            if err == 2:
                raise Exception('EK init failed', 'agrid incompatible with box size')
            elif err != 0:
                raise Exception('EK init failed', 'unknown error')




        # TODO:
        def checkpoint(self):
            raise Exception("Please implement this method in the pickle routine.")

        def add_reaction(self, shape):
            raise Exception("This method is not implemented yet.")

        def add_boundary(self, shape)
            raise Exception("This method is not implemented yet.")



# Suggested interface usage
# ek = Electrokinetics(params)
# system.actor.add(ek)

# pos_ions = ek.Species(charge=+1, concentration=0.06, D=0.001)
# neg_ions = ek.Species(charge=-1, concentration=0.07, D=0.001)

# ek.add_species(pos_ions, neg_ions)

# pos_ion[5,6,1].density = 0.08     #nodegrid 5,6,1


