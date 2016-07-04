include "myconfig.pxi"
cimport reaction
cimport globals
cimport utils
from highlander import ThereCanOnlyBeOne

IF CATALYTIC_REACTIONS:
    __reaction_is_initiated = False

    cdef class Reaction:

        def validate_params(self):
            if not isinstance(self._params["product_type"],int):
                raise ValueError("product_type has to be an int")
            if not isinstance(self._params["reactant_type"],int):
                raise ValueError("reactant_type has to be an int")
            if not isinstance(self._params["catalyzer_type"],int):
                raise ValueError("catalyzer_type has to be an int")
            if not isinstance(self._params["ct_range"],float):
                raise ValueError("ct_range has to be a float")
            if not isinstance(self._params["ct_rate"],float):
                raise ValueError("ct_rate has to be a float")
            if not isinstance(self._params["eq_rate"],float):
                raise ValueError("eq_rate has to be a float")
            if not isinstance(self._params["react_once"],bool):
                raise ValueError("react_once has to be a bool")
            if not isinstance(self._params["swap"],bool):
                raise ValueError("swap has to be a bool")

            if self._params["product_type"] < 0:
                raise ValueError("product_type has to be positive!")
            if self._params["reactant_type"] < 0:
                raise ValueError("reactant_type has to be positive!")
            if self._params["catalyzer_type"] < 0:
                raise ValueError("catalyzer_type has to be positive!")
            if( self._params["product_type"] == self._params["reactant_type"] or
                self._params["product_type"] == self._params["catalyzer_type"] or
                self._params["catalyzer_type"] == self._params["reactant_type"] ):
                raise Exception("One particle type cannot be part of more than one reaction species!")
            if self._params["ct_range"] < 0:
                raise ValueError("Negative equilibrium reaction rate contstant is not allowed!")
            if self._params["ct_rate"] < 0:
                raise ValueError("Negative catalytic reaction rate constant is not allowed!")
            if (self._params["eq_rate"] < 0 and abs(self._params["eq_rate"] + 1.0) > 0.001 ):
                raise ValueError("Negative equilibrium reaction rate contstant is not allowed!")

        def valid_keys(self):
            return "product_type", "reactant_type", "catalyzer_type", "ct_range", "ct_rate", "eq_rate", "react_once", "swap"

        def required_keys(self):
            return "product_type", "reactant_type", "catalyzer_type", "ct_range", "ct_rate"

        def default_params(self):
            return {"product_type": None,
                    "reactant_type": None,
                    "catalyzer_type": None,
                    "ct_range": None,
                    "ct_rate": None,
                    "eq_rate": 0.0,
                    "react_once": False,
                    "swap": False}

        def _set_params_in_es_core(self):
            globals.reaction.product_type   = self._params["product_type"]
            globals.reaction.reactant_type  = self._params["reactant_type"]
            globals.reaction.catalyzer_type = self._params["catalyzer_type"]
            globals.reaction.range   = self._params["ct_range"]
            globals.reaction.ct_rate = self._params["ct_rate"]
            globals.reaction.eq_rate = self._params["eq_rate"]
            globals.reaction.sing_mult = int(self._params["react_once"])
            globals.reaction.swap = int(self._params["swap"])

        def _get_params_from_es_core(self):
            self._params["product_type"]   = globals.reaction.product_type
            self._params["reactant_type"]  = globals.reaction.reactant_type
            self._params["catalyzer_type"] = globals.reaction.catalyzer_type
            self._params["ct_range"] = globals.reaction.range
            self._params["ct_rate"]  = globals.reaction.ct_rate
            self._params["eq_rate"]  = globals.reaction.eq_rate
            self._params["react_once"] = bool(globals.reaction.sing_mult)
            self._params["swap"] = bool(globals.reaction.swap)

        def get_params(self):
            self._get_params_from_es_core()
            return self._params

        # __getstate__ and __setstate__ define the pickle interaction
        def __getstate__(self):
            self._get_params_from_es_core()
            odict = self._params.copy()
            return odict

        def __setstate__(self,params):
            self._params=params
            self._set_params_in_es_core()
            mpi_setup_reaction()

        def __init__(self, *args, **kwargs):
            """
            Initialize the reaction.  Keep in mind, that there may be
            only one reaction enabled.  There can be only one.
            """
            self._ct_rate = 0.0

            # There can only be one reaction
            global __reaction_is_initiated
            if __reaction_is_initiated == True:
                raise ThereCanOnlyBeOne(self.__class__)
            __reaction_is_initiated = True

            # Get default params
            self._params = self.default_params()

            # Check if required keys are given
            for k in self.required_keys():
                if k not in kwargs:
                    raise ValueError("At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__() + " got " + kwargs.__str__())

                self._params[k] = kwargs[k]

            # Get default params
            self.setup(*args, **kwargs)

        def setup(self, *args, **kwargs):
            """
            Collect the parameters and set them in the core.
            """

            # Check if parameters are complete
            for k in kwargs:
                if k in self.valid_keys():
                    self._params[k] = kwargs[k]
                else:
                    raise KeyError("%s is not a valid key" % k)

            # Check if parameters are valid
            self.validate_params()

            # Pass parameters to Core
            self._set_params_in_es_core()

            # Set up the reaction
            mpi_setup_reaction()

        def start(self):
            """
            Restart the reaction after it was stopped
            """
            if (self._ct_rate != 0.0):
                self._params["ct_rate"] = self._ct_rate
                self._ct_rate = 0.0
                self._set_params_in_es_core()
                mpi_setup_reaction()

        def stop(self):
            """
            Stop the reaction, i.e. set the reaction rate to 0.0
            """
            if (self._ct_rate == 0.0):
                self._ct_rate = self._params["ct_rate"]
                self._params["ct_rate"] = 0.0
                self._set_params_in_es_core()
                mpi_setup_reaction()

        def __str__(self):
            return str(self.get_params())
