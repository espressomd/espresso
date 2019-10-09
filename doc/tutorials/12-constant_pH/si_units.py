import numpy as np
# physical constants
import scipy.constants
#as const

class SIunits():
    def __init__(self, **kwargs):
        self.params = kwargs
        required_params = ['sigma_in_nm', 'T_in_K']
        for p in required_params:
            if(p not in self.params.keys()):
                raise ValueError("\n" + p
                                 + " is missing.\nParameters that need to be defined for SI units: "
                                 + str(required_params) + "\n")
        # conversion factor for the unit of length
        self.sigma_to_nm = self.params['sigma_in_nm']
        # conversion factor for the unit of concentration
        self.conc_N_sigma_to_mol_L = 1.0e-3 / \
            (scipy.constants.N_A * (self.params['sigma_in_nm'] * 1e-9)**3)
        pass

    def convert(self, from_value, from_unit, to_unit, exponent=1.0):
        if (from_unit == "sigma" and to_unit == "nm"):
            return from_value * self.sigma_to_nm
        if (from_unit == "nm" and to_unit == "sigma"):
            return from_value / self.sigma_to_nm
        if (from_unit == "mol/L" and to_unit == "N/sigma^3"):
            return from_value * self.conc_N_sigma_to_mol_L
        if (to_unit == "mol/L" and from_unit == "N/sigma^3"):
            return from_value / self.conc_N_sigma_to_mol_L
        raise ValueError(
            "conversion from {} to {} is not implemented\n".format(
                from_unit, to_unit))
