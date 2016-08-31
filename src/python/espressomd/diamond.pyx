from __future__ import print_function, absolute_import
include "myconfig.pxi"
cdef class Diamond:
    def __init__(self, *args, **kwargs):
        self._params = self.default_params()
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError("At least the following keys have to be given as keyword arguments: " +self.required_keys().__str__() + " got " + kwargs.__str__())
        for k in kwargs:
            if k in self.valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a vaild key" % k)

        self.validate_params()
        self._set_params_in_es_core()

    def default_params(self):
        return {"a": 0.0, "bond_length": 0.0, "MPC": 0, "N_CI": 0, "val_nodes": 0.0, "val_cM": 0.0, "val_CI": 0.0, "cM_dist": 1, "nonet": 0}

    def required_keys(self):
        return "a", "bond_length", "MPC"

    def valid_keys(self):
        return "a", "bond_length", "MPC", "N_CI", "val_nodes", "val_cM", "val_CI", "cM_dist", "nonet"

    def validate_params(self):
        valid_keys = self.valid_keys()
        if(self._params[valid_keys[1]] < 0):
            raise ValueError("Bond length must be positive got: ",
                             self._params[valid_keys[1]])
        if(self._params[valid_keys[2]] < 0 or not isinstance(self._params[valid_keys[2]], int)):
            raise ValueError(
                "Monomers per chain must be positive integer got: ", self._params[valid_keys[2]])
        if(self._params[valid_keys[3]] < 0 or not isinstance(self._params[valid_keys[3]], int)):
            raise ValueError(
                "The number of counterions must be integer got:", self._params[valid_keys[3]])
        if(not isinstance(self._params[valid_keys[4]], float)):
            raise ValueError(
                "The charge of the nodes must be double got", self._params[valid_keys[4]])
        if(not isinstance(self._params[valid_keys[5]], float)):
            raise ValueError(
                "The charge of the monomers must be double got", self._params[valid_keys[5]])
        if(not isinstance(self._params[valid_keys[6]], float)):
            raise ValueError(
                "The charge of the counterions must be double got", self._params[valid_keys[6]])
        if(not isinstance(self._params[valid_keys[7]], int)):
            raise ValueError(
                "The distance between two charged monomers' indices must be integer ", self._params[valid_keys[7]])
        if(self._params[valid_keys[7]] == "nonet"):
            self._params[valid_keys[7]] = 1
        if(n_bonded_ia == 0):
            raise ValueError(
                "Please define a bonded interaction [0] before setting up polymers!")

    def __set_params_in_es_core(self):
        return diamondC(self._params["a"], self._params["bond_length"], self._params["MPC"], self._params["N_CI"], self._params["val_nodes"], self._params["val_cM"], self._params["val_CI"], self._params["cM_dist"], self._params["nonet"])

    def _set_params_in_es_core(self):
        tmp_try = self.__set_params_in_es_core()
        if(tmp_try == -3):
            raise Exception(
                "Failed upon creating one of the monomers in Espresso!\nAborting...\n")
        elif(tmp_try >= 0):
            print(tmp_try)
        else:
            raise Exception("Unknown Error")
