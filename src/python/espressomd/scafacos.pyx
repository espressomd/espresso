from __future__ import print_function, absolute_import
from espressomd.actors cimport *
from libcpp.string cimport string  # import std::string
cimport electrostatics
cimport magnetostatics

from espressomd.utils cimport handle_errors

include "myconfig.pxi"

# Interface to the scafacos library. These are the methods shared between
# dipolar and electrostatics methods

IF SCAFACOS == 1:
    class ScafacosConnector(Actor):

        def valid_keys(self):
            return "method_name", "method_params", "bjerrum_length"

        def required_keys(self):
            return "method_name", "method_params", "bjerrum_length"

        def validate_params(self):
            return True

        def _get_params_from_es_core(self):
            # Parameters are returned as strings
            # First word is method name, rest are key value pairs
            # which we convert to a dict
            p = get_method_and_parameters().split(" ")
            res = {}
            res["method_name"] = p[0]
            for i in range((len(p) - 1) / 2):
                res[p[2 * i + 1]] = p[2 * i + 2]
            return res

        def _set_params_in_es_core(self):
            # Verify that scafacos is not used for elecrostatics and dipoles
            # at the same time
            IF ELECTROSTATICS == 1:
                if self.dipolar and <int>electrostatics.coulomb.method ==<int>electrostatics.COULOMB_SCAFACOS:
                    raise Exception("Scafacos cannot be used for dipoles and charges at the same time")

            IF DIPOLES == 1:
                if not self.dipolar and <int>magnetostatics.coulomb.Dmethod ==<int>magnetostatics.DIPOLAR_SCAFACOS:
                    raise Exception("Scafacos cannot be used for dipoles and charges at the same time")


            # Convert the parameter dictionary to a list of strings
            method_params = self._params["method_params"]
            param_string = ""
            for k in method_params:
                param_string += k + " " + str(method_params[k]) + " "
            # Remove trailing whitespace
            param_string = param_string[0:-1]
            param_string = param_string.replace(" ", ",")

            set_parameters(self._params["method_name"],
                           param_string, self.dipolar)
            handle_errors("Scafacos not initialized.")

        def default_params(self):
            return {}
