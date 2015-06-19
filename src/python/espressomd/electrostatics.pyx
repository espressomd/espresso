#
# Copyright (C) 2013,2014 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
include "myconfig.pxi"
import numpy as np
from actors import Actor

class ElectrostaticInteraction(Actor):

    def _tune(self):
        raise Exception(
            "Subclasses of ElectrostaticInteraction must define the _tune() method or chosen method does not support tuning.")



IF COULOMB_DEBYE_HUECKEL:
    class CDH(ElectrostaticInteraction):
        def validateParams(self):
            if (self._params["bjerrum_length"] <= 0):
                raise ValueError("Bjerrum_length should be a positive double")
            if (self._params["kappa"] < 0):
                raise ValueError("kappa should be a non-negative double")
            if (self._params["r_cut"] < 0):
                raise ValueError("r_cut should be a non-negative double")
            if (self._params["eps_int"] <= 0):
                raise ValueError("eps_int should be a positive double")
            if (self._params["eps_ext"] <= 0):
                raise ValueError("eps_ext should be a positive double")
            if (self._params["r0"] <= 0 or self._params["r0"] >= self._params["r1"]):
                raise ValueError("r0 should be a positive double smaller than r1")
            if (self._params["r1"] <= 0 or self._params["r1"] >= self._params["r_cut"]):
                raise ValueError("r1 should be a positive double larger than r0 and smaller than r_cut")
            if (self._params["alpha"] < 0):
                raise ValueError("alpha should be a non-negative double")
        
        def validKeys(self):
            return "bjerrum_length","kappa","r_cut","eps_int","eps_ext","r0","r1","alpha"
        
        def requiredKeys(self):
            return "bjerrum_length","kappa","r_cut","eps_int","eps_ext","r0","r1","alpha"
        
        def _setParamsInEsCore(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            dh_set_params_cdh(self._params["kappa"], self._params["r_cut"], self._params["eps_int"], self._params["eps_ext"], self._params["r0"], self._params["r1"], self._params["alpha"])
        
        def _getParamsFromEsCore(self):
            params = {}
            params.update(dh_params)
            return params
        
        def _activateMethod(self):
            coulomb.method = COULOMB_DH

        def defaultParams(self):
            return {"bjerrum_length": -1,
                    "kappa": -1,
                    "r_cut": -1,
                    "eps_int": -1,
                    "eps_ext": -1,
                    "r0": -1,
                    "r1": -1,
                    "alpha": -1}

ELSE:
    IF ELECTROSTATICS:
        class DH(ElectrostaticInteraction):
            def validateParams(self):
                if (self._params["bjerrum_length"] <= 0):
                    raise ValueError("Bjerrum_length should be a positive double")
                if (self._params["kappa"] <0):
                    raise ValueError("kappa should be a non-negative double")
                if (self._params["r_cut"] <0):
                    raise ValueError("r_cut should be a non-negative double")
            
            def validKeys(self):
                return "bjerrum_length","kappa","r_cut"
        
            def requiredKeys(self):   
                return "bjerrum_length","kappa","r_cut"
        
            def _setParamsInEsCore(self):
                coulomb_set_bjerrum(self._params["bjerrum_length"])
                dh_set_params(self._params["kappa"], self._params["r_cut"])
            
            def _getParamsFromEsCore(self):
                params = {}
                params.update(dh_params)
                return params
            
            def _activateMethod(self):
                coulomb.method = COULOMB_DH
                self._setParamsInEsCore()

            def defaultParams(self):
                return {"bjerrum_length": -1,
                        "kappa": -1,
                        "r_cut": -1}


IF P3M == 1:
    class P3M(ElectrostaticInteraction):

        def validateParams(self):
            default_params = self.defaultParams()
            if not (self._params["bjerrum_length"] > 0.0):
                raise ValueError("Bjerrum_length should be a positive double")

            if not (self._params["r_cut"] >= 0 or self._params["r_cut"] == default_params["r_cut"]):
                raise ValueError("P3M r_cut has to be >=0")

            if not (isinstance(self._params["mesh"], int) or len(self._params["mesh"])):
                raise ValueError(
                    "P3M mesh has to be an integer or integer list of length 3")

            if (isinstance(self._params["mesh"], basestring) and len(self._params["mesh"]) == 3):
                if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                    raise ValueError(
                        "P3M requires an even number of mesh points in all directions")

            if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                raise ValueError(
                    "P3M cao has to be an integer between -1 and 7")

            if not (self._params["accuracy"] > 0):
                raise ValueError("P3M accuracy has to be positive")

            if self._params["epsilon"] == "metallic":
                self._params = 0.0

            if not (isinstance(self._params["epsilon"], float) or self._params["epsilon"] == "metallic"):
                raise ValueError("epsilon should be a double or 'metallic'")

            if not (self._params["inter"] == default_params["inter"] or self._params["inter"] > 0):
                raise ValueError("inter should be a positive integer")

            if not (self._params["mesh_off"] == default_params["mesh_off"] or len(self._params) != 3):
                raise ValueError(
                    "mesh_off should be a list of length 3 and values between 0.0 and 1.0")

        def validKeys(self):
            return "alpha_L", "r_cut_iL", "mesh", "mesh_off", "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai", "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "bjerrum_length", "tune"

        def requiredKeys(self):
            return ["bjerrum_length", "accuracy"]

        def defaultParams(self):
            return {"cao": -1,
                    "inter": -1,
                    "r_cut": -1,
                    "accuracy": -1,
                    "mesh": [-1, -1, -1],
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True}

        def _getParamsFromEsCore(self):
            params = {}
            params.update(p3m.params)
            params["bjerrum_length"] = coulomb.bjerrum
            params["tune"] = self._params["tune"]
            return params

        def _setParamsInEsCore(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            p3m_set_ninterpol(self._params["inter"])
            python_p3m_set_mesh_offset(self._params["mesh_off"])
            python_p3m_set_params(self._params["r_cut"], self._params["mesh"], self._params[
                                  "cao"], self._params["alpha"], self._params["accuracy"])
            p3m_set_eps(self._params["epsilon"])

        def _tune(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            p3m_set_eps(self._params["epsilon"])
            python_p3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params[
                                       "cao"], -1.0, self._params["accuracy"], self._params["inter"])
            resp, log = python_p3m_adaptive_tune()
            if resp:
                raise Exception(
                    "failed to tune P3M parameters to required accuracy")
            print log
            self._params.update(self._getParamsFromEsCore())

        def _activateMethod(self):
            if self._params["tune"]:
                self._tune()

            self._setParamsInEsCore()

    IF CUDA:
        class P3M_GPU(ElectrostaticInteraction):

            def validateParams(self):
                default_params = self.defaultParams()
                if not (self._params["bjerrum_length"] > 0.0):
                    raise ValueError(
                        "Bjerrum_length should be a positive double")

                if not (self._params["r_cut"] >= 0 or self._params["r_cut"] == default_params["r_cut"]):
                    raise ValueError("P3M r_cut has to be >=0")

                if not (isinstance(self._params["mesh"], int) or len(self._params["mesh"])):
                    raise ValueError(
                        "P3M mesh has to be an integer or integer list of length 3")

                if (isinstance(self._params["mesh"], basestring) and len(self._params["mesh"]) == 3):
                    if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                        raise ValueError(
                            "P3M requires an even number of mesh points in all directions")

                if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                    raise ValueError(
                        "P3M cao has to be an integer between -1 and 7")

                if not (self._params["accuracy"] > 0):
                    raise ValueError("P3M accuracy has to be positive")

                # if self._params["epsilon"] == "metallic":
                #  self._params = 0.0

                if not (isinstance(self._params["epsilon"], float) or self._params["epsilon"] == "metallic"):
                    raise ValueError(
                        "epsilon should be a double or 'metallic'")

                if not (self._params["inter"] == default_params["inter"] or self._params["inter"] > 0):
                    raise ValueError("inter should be a positive integer")

                if not (self._params["mesh_off"] == default_params["mesh_off"] or len(self._params) != 3):
                    raise ValueError(
                        "mesh_off should be a list of length 3 and values between 0.0 and 1.0")

            def validKeys(self):
                return "alpha_L", "r_cut_iL", "mesh", "mesh_off", "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai", "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "bjerrum_length", "tune"

            def requiredKeys(self):
                return ["bjerrum_length", "accuracy"]

            def defaultParams(self):
                return {"cao": -1,
                        "inter": -1,
                        "r_cut": -1,
                        "accuracy": -1,
                        "mesh": [-1, -1, -1],
                        "epsilon": 0.0,
                        "mesh_off": [-1, -1, -1],
                        "tune": True,
                        "box": [-1, -1, -1],
                        "alpha": -1}

            def _getParamsFromEsCore(self):
                params = {}
                params.update(p3m.params)
                params["bjerrum_length"] = coulomb.bjerrum
                params["tune"] = self._params["tune"]
                params["box"] = self.system.box_l
                return params

            def _tune(self):
                python_p3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params[
                                           "cao"], -1.0, self._params["accuracy"], self._params["inter"])
                resp, log = python_p3m_adaptive_tune()
                if resp:
                    raise Exception(
                        "failed to tune P3M parameters to required accuracy")
                print log
                self._params.update(self._getParamsFromEsCore())

            def _activateMethod(self):
                coulomb.method = COULOMB_P3M_GPU
                if self._params["tune"]:
                    coulomb_set_bjerrum(self._params["bjerrum_length"])
                    self._tune()

                coulomb_set_bjerrum(self._params["bjerrum_length"])
                self._setParamsInEsCore()
                python_p3m_gpu_init(self._params)

            def _setParamsInEsCore(self):
                python_p3m_set_params(self._params["r_cut"], self._params["mesh"], self._params[
                                      "cao"], self._params["alpha"], self._params["accuracy"])
                p3m_set_eps(self._params["epsilon"])
                coulomb_set_bjerrum(self._params["bjerrum_length"])
                p3m_set_ninterpol(self._params["inter"])
                python_p3m_set_mesh_offset(self._params["mesh_off"])
