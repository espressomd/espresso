#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from __future__ import print_function, absolute_import
from cython.operator cimport dereference
include "myconfig.pxi"
from espressomd cimport actors
from . import actors
cimport globals
import numpy as np
from .scafacos import *
from . cimport scafacos


cdef class ElectrostaticInteraction(actors.Actor):
    def _tune(self):
        raise Exception(
            "Subclasses of ElectrostaticInteraction must define the _tune() method or chosen method does not support tuning.")

    def _set_params_in_es_core(self):
        raise Exception(
            "Subclasses of ElectrostaticInteraction must define the _set_params_in_es_core() method.")

    def _deactivate_method(self):
        coulomb.method = COULOMB_NONE
        mpi_bcast_coulomb_params()

    def Tune(self,subsetTuneParams=None):
        
        #Override default parmas with subset given by user
        tuneParams=self.default_params()
        if not subsetTuneParams==None:        
            for k in subsetTuneParams.iterkeys():
                if k not in self.valid_keys():
                    raise ValueError(k + " is not a valid parameter")
            tuneParams.update(subsetTuneParams)

        #If param is 'required', it was set before, so don't change it
        #Do change it if it's given to Tune() by user
        for param in tuneParams.iterkeys():
            if not param in self.required_keys() or (not subsetTuneParams==None and param in subsetTuneParams.keys()):
                self._params[param]=tuneParams[param]
        self._set_params_in_es_core()
        self._tune()

IF COULOMB_DEBYE_HUECKEL:
    cdef class CDH(ElectrostaticInteraction):
        def validate_params(self):
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
                raise ValueError(
                    "r0 should be a positive double smaller than r1")
            if (self._params["r1"] <= 0 or self._params["r1"] >= self._params["r_cut"]):
                raise ValueError(
                    "r1 should be a positive double larger than r0 and smaller than r_cut")
            if (self._params["alpha"] < 0):
                raise ValueError("alpha should be a non-negative double")

        def valid_keys(self):
            return "bjerrum_length", "kappa", "r_cut", "eps_int", "eps_ext", "r0", "r1", "alpha"

        def required_keys(self):
            return "bjerrum_length", "kappa", "r_cut", "eps_int", "eps_ext", "r0", "r1", "alpha"

        def _set_params_in_es_core(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            dh_set_params_cdh(self._params["kappa"], self._params["r_cut"], self._params[
                              "eps_int"], self._params["eps_ext"], self._params["r0"], self._params["r1"], self._params["alpha"])

        def _get_params_from_es_core(self):
            params = {}
            params.update(dh_params)
            return params

        def _activate_method(self):
            coulomb.method = COULOMB_DH

        def default_params(self):
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
        cdef class DH(ElectrostaticInteraction):
            def validate_params(self):
                if (self._params["bjerrum_length"] <= 0):
                    raise ValueError(
                        "Bjerrum_length should be a positive double")
                if (self._params["kappa"] < 0):
                    raise ValueError("kappa should be a non-negative double")
                if (self._params["r_cut"] < 0):
                    raise ValueError("r_cut should be a non-negative double")

            def valid_keys(self):
                return "bjerrum_length", "kappa", "r_cut"

            def required_keys(self):
                return "bjerrum_length", "kappa", "r_cut"

            def _set_params_in_es_core(self):
                coulomb_set_bjerrum(self._params["bjerrum_length"])
                dh_set_params(self._params["kappa"], self._params["r_cut"])

            def _get_params_from_es_core(self):
                params = {}
                params.update(dh_params)
                return params

            def _activate_method(self):
                coulomb.method = COULOMB_DH
                self._set_params_in_es_core()

            def default_params(self):
                return {"bjerrum_length": -1,
                        "kappa": -1,
                        "r_cut": -1}


IF P3M == 1:
    cdef class P3M(ElectrostaticInteraction):

        def validate_params(self):
            default_params = self.default_params()
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

            if not (self._params["accuracy"] >= 0):
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

            if not (self._params["alpha"] == default_params["alpha"] or self._params["alpha"] > 0):
                raise ValueError(
                    "alpha should be positive")

        def valid_keys(self):
            return "alpha_L", "r_cut_iL", "mesh", "mesh_off", "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai", "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "bjerrum_length", "tune"

        def required_keys(self):
            return ["bjerrum_length"]

        def default_params(self):
            return {"cao": 0,
                    "inter": -1,
                    "r_cut": -1,
                    "alpha": 0,
                    "accuracy": 0,
                    "mesh": [0, 0, 0],
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(p3m.params)
            params["bjerrum_length"] = coulomb.bjerrum
            params["tune"] = self._params["tune"]
            return params

        def _set_params_in_es_core(self):
            #Sets lb, bcast, resets vars to zero if lb=0
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            #Sets cdef vars and calls p3m_set_params() in core 
            python_p3m_set_params(self._params["r_cut"],
                        self._params["mesh"], self._params["cao"],
                        self._params["alpha"], self._params["accuracy"])
            #p3m_set_params()  -> set r_cuts, mesh, cao, validates sanity, bcasts 
            #Careful: bcast calls on_coulomb_change(), which calls p3m_init(),
            #         which resets r_cut if lb is zero. OK.
            #Sets eps, bcast
            p3m_set_eps(self._params["epsilon"])
            #Sets ninterpol, bcast
            p3m_set_ninterpol(self._params["inter"])
            python_p3m_set_mesh_offset(self._params["mesh_off"])

        def _tune(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            python_p3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params[
                                       "cao"], -1.0, self._params["accuracy"], self._params["inter"])
            resp = python_p3m_adaptive_tune()
            if resp:
                raise Exception(
                    "failed to tune P3M parameters to required accuracy")
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()

    IF CUDA:
        cdef class P3M_GPU(ElectrostaticInteraction):

            def validate_params(self):
                default_params = self.default_params()
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

                if not (self._params["accuracy"] >= 0):
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

            def valid_keys(self):
                return "alpha_L", "r_cut_iL", "mesh", "mesh_off", "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai", "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "bjerrum_length", "tune"

            def required_keys(self):
                return ["bjerrum_length"]

            def default_params(self):
                return {"cao": 0,
                        "inter": -1,
                        "r_cut": -1,
                        "alpha": 0,
                        "accuracy": 0,
                        "mesh": [0, 0, 0],
                        "epsilon": 0.0,
                        "mesh_off": [-1, -1, -1],
                        "tune": True}

            def _get_params_from_es_core(self):
                params = {}
                params.update(p3m.params)
                params["bjerrum_length"] = coulomb.bjerrum
                params["tune"] = self._params["tune"]
                params["box"] = self.system.box_l
                return params

            def _tune(self):
                python_p3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params[
                                           "cao"], -1.0, self._params["accuracy"], self._params["inter"])
                resp = python_p3m_adaptive_tune()
                if resp:
                    raise Exception(
                        "failed to tune P3M parameters to required accuracy")
                self._params.update(self._get_params_from_es_core())

            def _activate_method(self):
                self._set_params_in_es_core()
                coulomb.method = COULOMB_P3M_GPU
                #python_p3m_gpu_init(self._params)
                if self._params["tune"]:
                    self._tune()

                self._set_params_in_es_core()

            def _set_params_in_es_core(self):
                python_p3m_set_params(self._params["r_cut"], self._params["mesh"], self._params[
                                      "cao"], self._params["alpha"], self._params["accuracy"])
                p3m_set_eps(self._params["epsilon"])
                coulomb_set_bjerrum(self._params["bjerrum_length"])
                p3m_set_ninterpol(self._params["inter"])
                python_p3m_set_mesh_offset(self._params["mesh_off"])


IF ELECTROSTATICS and CUDA and EWALD_GPU:
    cdef class EwaldGpu(ElectrostaticInteraction):
        cdef EwaldgpuForce * thisptr
        cdef EspressoSystemInterface * interface
        cdef char * log
        cdef int resp

        def __cinit__(self):
            self.interface = EspressoSystemInterface._Instance()
            default_params = self.default_params()
            self.thisptr = new EwaldgpuForce(dereference(self.interface), default_params["rcut"], default_params["num_kx"], default_params["num_ky"], default_params["num_kz"], default_params["alpha"])

        def __dealloc__(self):
            del self.thisptr

        def valid_keys(self):
            return "bjerrum_length", "rcut", "num_kx", "num_ky", "num_kz", "K_max", "alpha", "accuracy", "precision", "time_calc_steps"

        def default_params(self):
            return {"bjerrum_length": -1,
                    "rcut": -1,
                    "num_kx": -1,
                    "num_ky": -1,
                    "num_kz": -1,
                    "alpha": -1,
                    "accuracy": 0,
                    "precision": -1,
                    "isTuned": False,
                    "isTunedFlag": False,
                    "K_max": -1,
                    "time_calc_steps": -1}

        def validate_params(self):
            default_params = self.default_params()
            if self._params["bjerrum_length"] <= 0.0 and self._params["bjerrum_length"] != default_params["bjerrum_length"]:
                raise ValueError("Bjerrum_length should be a positive double")
            if isinstance(self._params["K_max"], (list, np.ndarray)):
                if isinstance(self._params["K_max"], int) and len(self._params["K_max"]) == 3:
                    if self._params["K_max"][0] < 0 or self._params["K_max"][1] < 0 or self._params["K_max"][2] < 0:
                        raise ValueError(
                            "K_max has to be a positive integer or a list of three positive integers")
            elif self._params["K_max"] < 0:
                raise ValueError(
                    "K_max has to be a positive integer or a list of three positive integers")
            if self._params["rcut"] < 0 and self._params["rcut"] != default_params["rcut"]:
                raise ValueError("rcut should be a positive float")
            if self._params["accuracy"] < 0 and self._params["accuracy"] != default_params["accuracy"]:
                raise ValueError("accuracy has to be a positive double")
            if self._params["precision"] < 0 and self._params["precision"] != default_params["precision"]:
                raise ValueError("precision has to be a positive double")
            if self._params["alpha"] < 0 and self._params["alpha"] != default_params["alpha"]:
                raise ValueError("alpha has to be a positive double")

        def required_keys(self):
            return "bjerrum_length", "accuracy", "precision", "K_max"

        def _tune(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            default_params = self.default_params()
            if self._params["time_calc_steps"] == default_params["time_calc_steps"]:
                self._params[
                    "time_calc_steps"] = self.thisptr.determine_calc_time_steps()
            if isinstance(self._params["K_max"], (list, np.ndarray)):
                self.thisptr.set_params(self._params["rcut"], self._params["K_max"][0], self._params[
                                        "K_max"][1], self._params["K_max"][2], self._params["alpha"])
            else:
                self.thisptr.set_params_tune(self._params["accuracy"], self._params[
                    "precision"], self._params["K_max"], self._params["time_calc_steps"])
            resp = self.thisptr.adaptive_tune( & self.log, dereference(self.interface))
            if resp != 0:
                print(self.log)

        def _set_params_in_es_core(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            self.thisptr.set_params(self._params["rcut"], self._params[
                                    "num_kx"], self._params["num_ky"], self._params["num_kz"], self._params["alpha"])

        def _get_params_from_es_core(self):
            params = {}
            params.update(ewaldgpu_params)
            params["bjerrum_length"] = coulomb.bjerrum
            return params

        def _activate_method(self):
            self._set_params_in_es_core()
            coulomb.method = COULOMB_EWALD_GPU
            if not self._params["isTuned"]:
                self._tune()
                self._params["isTuned"] = True

            self._set_params_in_es_core()

IF ELECTROSTATICS:
    cdef class MMM1D(electrostatics.ElectrostaticInteraction):

        def validate_params(self):
            default_params = self.default_params()
            if self._params["bjerrum_length"] < 0:
                raise ValueError("Bjerrum_length should be a positive double")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["far_switch_radius_2"] < 0 and self._params["far_switch_radius_2"] != default_params["far_switch_radius_2"]:
                raise ValueError("switch radius shoulb be a positive double")
            if self._params["far_switch_radius"] < 0 and self._params["far_switch_radius"] != default_params["far_switch_radius"]:
                raise ValueError("switch radius shoulb be a positive double")
            if self._params["bessel_cutoff"] < 0 and self._params["bessel_cutoff"] != default_params["bessel_cutoff"]:
                raise ValueError("bessel_cutoff should be a positive integer")

        def default_params(self):
            return {"bjerrum_length": -1,
                    "maxPWerror": -1,
                    "far_switch_radius_2": -1,
                    "far_switch_radius": -1,
                    "bessel_cutoff": -1,
                    "tune": True}

        def valid_keys(self):
            return "bjerrum_length", "maxPWerror", "far_switch_radius", "bessel_cutoff", "tune"

        def required_keys(self):
            return ["bjerrum_length", "maxPWerror"]

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm1d_params)
            params["far_switch_radius"] = np.sqrt(
                params["far_switch_radius_2"])
            params["bjerrum_length"] = coulomb.bjerrum
            return params

        def _set_params_in_es_core(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            if self._params["far_switch_radius"] == -1:
                self._params["far_switch_radius_2"] = -1
            else:
                self._params["far_switch_radius_2"] = self._params[
                    "far_switch_radius"] * self._params["far_switch_radius"]
            MMM1D_set_params(
                self._params["far_switch_radius_2"], self._params["maxPWerror"])

        def _tune(self):
            cdef int resp
            resp = pyMMM1D_tune()
            if resp:
                raise Exception("failed to tune mmm1d ")
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            coulomb.method = COULOMB_MMM1D
            self._set_params_in_es_core()
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()

IF ELECTROSTATICS and MMM1D_GPU:
    cdef class MMM1D_GPU(ElectrostaticInteraction):
        cdef Mmm1dgpuForce * thisptr
        cdef EspressoSystemInterface * interface
        cdef char * log
        cdef int resp

        def __cinit__(self):
            self.interface = EspressoSystemInterface._Instance()
            default_params = self.default_params()
            self.thisptr = new Mmm1dgpuForce(dereference(self.interface), 0.0, default_params["maxPWerror"])
            self.interface.update()
            self.interface.requestRGpu()

        def __dealloc__(self):
            del self.thisptr

        def validate_params(self):
            default_params = self.default_params()
            if self._params["bjerrum_length"] < 0:
                raise ValueError("Bjerrum_length should be a positive double")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["far_switch_radius"] < 0 and self._params["far_switch_radius"] != default_params["far_switch_radius"]:
                raise ValueError("switch radius shoulb be a positive double")
            if self._params["bessel_cutoff"] < 0 and self._params["bessel_cutoff"] != default_params["bessel_cutoff"]:
                raise ValueError("bessel_cutoff should be a positive integer")

        def default_params(self):
            return {"bjerrum_length": -1,
                    "maxPWerror": -1.0,
                    "far_switch_radius": -1.0,
                    "far_switch_radius_2": -1.0,
                    "bessel_cutoff": -1,
                    "tune": True}

        def valid_keys(self):
            return "bjerrum_length", "maxPWerror", "far_switch_radius", "bessel_cutoff", "tune"

        def required_keys(self):
            return ["bjerrum_length", "maxPWerror"]

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm1d_params)
            params["far_switch_radius"] = np.sqrt(
                params["far_switch_radius_2"])
            params["bjerrum_length"] = coulomb.bjerrum
            return params

        def _set_params_in_es_core(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            default_params = self.default_params()
            if self._params["far_switch_radius"] == default_params["far_switch_radius"]:
                self._params["far_switch_radius_2"] = -1
            else:
                self._params["far_switch_radius_2"] = self._params[
                    "far_switch_radius"] * self._params["far_switch_radius"]

            self.thisptr.set_params(globals.box_l[2], globals.temperature * coulomb.bjerrum, self._params[
                                    "maxPWerror"], self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _tune(self):
            self.thisptr.setup(dereference(self.interface))
            self.thisptr.tune(dereference(self.interface), self._params[
                              "maxPWerror"], self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _activate_method(self):
            self._set_params_in_es_core()
            coulomb.method = COULOMB_MMM1D_GPU
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()

IF ELECTROSTATICS:
    cdef class MMM2D(ElectrostaticInteraction):
        def validate_params(self):
            default_params = self.default_params()
            if self._params["bjerrum_length"] < 0:
                raise ValueError("Bjerrum_length should be a positive double")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["dielectric"] == 1 and (self._params["top"] < 0 or self._params["mid"] < 0 or self._params["bot"] < 0):
                raise ValueError("Dielectric constants should be > 0!")
            if self._params["dielectric_contrast_on"] == 1 and (self._params["delta_mid_top"] == default_params["delta_mid_top"] or self._params["delta_mid_bot"] == default_params["delta_mid_bot"]):
                raise ValueError("Dielectric constrast not set!")
            if self._params["capacitor"] == 1 and self._params["pot_diff"] == default_params["pot_diff"]:
                raise ValueError("Potential difference not set!")
            if self._params["dielectric"] == 1 and self._params["dielectric_contrast_on"] == 1:
                raise ValueError(
                    "dielectric and dielectric_contrast are mutually exclusive!")
            if self._params["dielectric"] == 1 and self._params["capacitor"] == 1:
                raise ValueError(
                    "dielectric and constant potential are mutually exclusive")
            if self._params["dielectric_contrast_on"] == 1 and self._params["capacitor"] == 1:
                raise ValueError(
                    "dielectric contrast and constant potential are mutually exclusive")

        def default_params(self):
            return {"bjerrum_length": -1,
                    "maxPWerror": -1,
                    "far_cut": -1,
                    "top": 0,
                    "mid": 0,
                    "bot": 0,
                    "dielectric": 0,
                    "top": 0,
                    "mid": 0,
                    "bot": 0,
                    "dielectric_contrast_on": 0,
                    "capacitor": 0,
                    "delta_mid_top": 0,
                    "delta_mid_bot": 0,
                    "pot_diff": 0}

        def required_keys(self):
            return ["bjerrum_length", "maxPWerror"]

        def valid_keys(self):
            return "bjerrum_length", "maxPWerror", "top", "mid", "bot", "delta_mid_top", "delta_mid_bot", "pot_diff", "dielectric", "dielectric_contrast_on", "capacitor", "far_cut"

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm2d_params)
            params["bjerrum_length"] = coulomb.bjerrum
            if params["dielectric_contrast_on"] or params["const_pot_on"]:
                params["dielectric"] = 0
            else:
                params["dielectric"] = 1
            return params

        def _set_params_in_es_core(self):
            coulomb_set_bjerrum(self._params["bjerrum_length"])
            if self._params["dielectric"]:
                self._params["delta_mid_top"] = (self._params[
                                                 "mid"] - self._params["top"]) / (self._params["mid"] + self._params["top"])
                self._params["delta_mid_bot"] = (self._params[
                                                 "mid"] - self._params["bot"]) / (self._params["mid"] + self._params["bot"])

            if self._params["capacitor"]:
                self._params["delta_mid_top"] = -1
                self._params["delta_mid_bot"] = -1
                self._params["const_pot_on"] = 1

            print(MMM2D_set_params(self._params["maxPWerror"], self._params["far_cut"], self._params["delta_mid_top"], self._params["delta_mid_bot"], self._params["capacitor"], self._params["pot_diff"]))

        def _activate_method(self):
            coulomb.method = COULOMB_MMM2D
            self._set_params_in_es_core()
            MMM2D_init()
            print(MMM2D_sanity_checks())

    IF SCAFACOS == 1:
        class Scafacos(ScafacosConnector, ElectrostaticInteraction):
            dipolar = False

            # Explicit constructor needed due to multiple inheritance
            def __init__(self, *args, **kwargs):
                actors.Actor.__init__(self, *args, **kwargs)

            def _activate_method(self):
                coulomb.method = COULOMB_SCAFACOS
                coulomb_set_bjerrum(self._params["bjerrum_length"])
                self._set_params_in_es_core()

            def default_params(self):
                return {}
