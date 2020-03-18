#
# Copyright (C) 2013-2019 The ESPResSo project
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
from cpython.exc cimport PyErr_CheckSignals, PyErr_SetInterrupt
include "myconfig.pxi"
from .utils cimport handle_errors, check_type_or_throw_except
from . cimport integrate

cdef class IntegratorHandle:
    """
    Provide access to all integrators.

    """

    cdef object _integrator

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        return self._integrator

    def __setstate__(self, state):
        self._integrator = state
        self._integrator._set_params_in_es_core()

    def get_state(self):
        """
        Return the integrator.

        """
        return self.__getstate__()

    def __init__(self):
        self.set_nvt()

    def __str__(self):
        return self.__class__.__name__ + \
            "(" + str(self._integrator.__class__.__name__) + ")"

    def run(self, *args, **kwargs):
        """
        Run the integrator.

        """
        return self._integrator.run(*args, **kwargs)

    def set_steepest_descent(self, *args, **kwargs):
        """
        Set the integration method to steepest descent
        (:class:`SteepestDescent`).

        .. seealso::
            :func:`espressomd.minimize_energy.steepest_descent`

        """
        self._integrator = SteepestDescent(*args, **kwargs)

    def set_vv(self):
        """
        Set the integration method to velocity Verlet, which is suitable for
        simulations in the NVT ensemble (:class:`VelocityVerlet`).

        """
        self._integrator = VelocityVerlet()

    def set_nvt(self):
        """
        Set the integration method to velocity Verlet, which is suitable for
        simulations in the NVT ensemble (:class:`VelocityVerlet`).

        """
        self._integrator = VelocityVerlet()

    def set_isotropic_npt(self, *args, **kwargs):
        """
        Set the integration method to a modified velocity Verlet designed for
        simulations in the NpT ensemble (:class:`VelocityVerletIsotropicNPT`).

        """
        self._integrator = VelocityVerletIsotropicNPT(*args, **kwargs)

    def set_brownian_dynamics(self):
        """
        Set the integration method to BD.

        """
        self._integrator = BrownianDynamics()


cdef class Integrator:
    """
    Integrator class.

    """
    cdef object _params

    def __getstate__(self):
        return self._params

    def __setstate__(self, params):
        self._params = params
        self._set_params_in_es_core()

    def _set_params_in_es_core(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of Integrator must define the _set_params_in_es_core() method.")

    def get_params(self):
        """Get integrator parameters.

        """
        return self.__getstate__()

    def __init__(self, *args, **kwargs):

        # Check if all required keys are given
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        self._params = self.default_params()
        self._params.update(kwargs)
        self.validate_params()
        self._set_params_in_es_core()

    def __str__(self):
        return self.__class__.__name__ + "(" + str(self.get_params()) + ")"

    def validate_params(self):
        """Check that parameters are valid.

        """
        return True

    def default_params(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of Integrator must define the default_params() method.")

    def valid_keys(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of Integrator must define the valid_keys() method.")

    def required_keys(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of Integrator must define the required_keys() method.")

    def run(self, steps=1, recalc_forces=False, reuse_forces=False):
        """
        Run the integrator.

        Parameters
        ----------
        steps : :obj:`int`
            Number of time steps to integrate.
        recalc_forces : :obj:`bool`, optional
            Recalculate the forces regardless of whether they are reusable.
        reuse_forces : :obj:`bool`, optional
            Reuse the forces from previous time step.

        """
        check_type_or_throw_except(steps, 1, int, "steps must be an int")
        assert steps >= 0, "steps has to be positive"
        check_type_or_throw_except(
            recalc_forces, 1, bool, "recalc_forces has to be a bool")
        check_type_or_throw_except(
            reuse_forces, 1, bool, "reuse_forces has to be a bool")

        _integrate(steps, recalc_forces, reuse_forces)

        if integrate.set_py_interrupt:
            PyErr_SetInterrupt()
            integrate.set_py_interrupt = False
            PyErr_CheckSignals()

        handle_errors("Encountered errors during integrate")


cdef class SteepestDescent(Integrator):
    """
    Steepest descent algorithm for energy minimization.

    Particles located at :math:`\\vec{r}_i` at integration step :math:`i` and
    experiencing a potential :math:`\mathcal{H}(\\vec{r}_i)` are displaced
    according to the equation:

    :math:`\\vec{r}_{i+1} = \\vec{r}_i - \\gamma\\nabla\mathcal{H}(\\vec{r}_i)`

    Parameters
    ----------
    f_max : :obj:`float`
        Convergence criterion. Minimization stops when the maximal force on
        particles in the system is lower than this threshold. Set this to 0
        when running minimization in a loop that stops when a custom
        convergence criterion is met.
    gamma : :obj:`float`
        Dampening constant.
    max_displacement : :obj:`float`
        Maximal allowed displacement per step. Typical values for a LJ liquid
        are in the range of 0.1% to 10% of the particle sigma.
    max_steps : :obj:`int`
        Maximal number of iterations. Gets updated at each call to :meth:`run`.

    """

    def default_params(self):
        return {"max_steps": 0}

    def valid_keys(self):
        """All parameters that can be set.

        """
        return {"f_max", "gamma", "max_displacement", "max_steps"}

    def required_keys(self):
        """Parameters that have to be set.

        """
        return {"f_max", "gamma", "max_displacement"}

    def validate_params(self):
        check_type_or_throw_except(
            self._params["f_max"], 1, float, "f_max must be a float")
        check_type_or_throw_except(
            self._params["gamma"], 1, float, "gamma must be a float")
        check_type_or_throw_except(
            self._params["max_displacement"], 1, float, "max_displacement must be a float")
        check_type_or_throw_except(
            self._params["max_steps"], 1, int, "max_steps must be an int")

    def _set_params_in_es_core(self):
        if integrate_set_steepest_descent(self._params["f_max"],
                                          self._params["gamma"],
                                          self._params["max_steps"],
                                          self._params["max_displacement"]):
            handle_errors(
                "Encountered errors setting up the steepest descent integrator")

    def run(self, steps=1, **kwargs):
        """
        Run the steepest descent.

        Parameters
        ----------
        steps : :obj:`int`
            Maximal number of time steps to integrate.

        Returns
        -------
        :obj:`int`
            Number of integrated steps.

        """
        check_type_or_throw_except(steps, 1, int, "steps must be an int")
        assert steps >= 0, "steps has to be positive"
        self._params["max_steps"] = steps

        integrated = mpi_steepest_descent(steps)

        handle_errors("Encountered errors during integrate")

        return integrated


cdef class VelocityVerlet(Integrator):
    """
    Velocity Verlet integrator, suitable for simulations in the NVT ensemble.

    """

    def default_params(self):
        return {}

    def valid_keys(self):
        """All parameters that can be set.

        """
        return {}

    def required_keys(self):
        """Parameters that have to be set.

        """
        return {}

    def validate_params(self):
        return True

    def _set_params_in_es_core(self):
        integrate_set_nvt()


IF NPT:
    cdef class VelocityVerletIsotropicNPT(Integrator):
        """
        Modified velocity Verlet integrator, suitable for simulations in the
        NpT ensemble with isotropic rescaling. Requires the NpT thermostat,
        activated with :meth:`espressomd.thermostat.Thermostat.set_npt`.

        Parameters
        ----------
        ext_pressure : :obj:`float`
            The external pressure.
        piston : :obj:`float`
            The mass of the applied piston.
        direction : (3,) array_like of :obj:`bool`, optional
            Select which dimensions are allowed to fluctuate by assigning
            them to ``True``.
        cubic_box : :obj:`bool`, optional
            If ``True``, a cubic box is assumed and the value of ``direction``
            will be ignored when rescaling the box. This is required e.g. for
            electrostatics and magnetostatics.
        """

        def default_params(self):
            return {"direction": (True, True, True), "cubic_box": False}

        def valid_keys(self):
            """All parameters that can be set.

            """
            return {"ext_pressure", "piston", "direction", "cubic_box"}

        def required_keys(self):
            """Parameters that have to be set.

            """
            return {"ext_pressure", "piston"}

        def validate_params(self):
            check_type_or_throw_except(
                self._params["ext_pressure"], 1, float, "ext_pressure must be a float")
            check_type_or_throw_except(
                self._params["piston"], 1, float, "piston must be a float")
            check_type_or_throw_except(
                self._params["direction"], 3, int, "direction must be an array-like of 3 bools")
            check_type_or_throw_except(
                self._params["cubic_box"], 1, int, "cubic_box must be a bool")

        def _set_params_in_es_core(self):
            if integrate_set_npt_isotropic(self._params["ext_pressure"],
                                           self._params["piston"],
                                           self._params["direction"][0],
                                           self._params["direction"][1],
                                           self._params["direction"][2],
                                           self._params["cubic_box"]):
                handle_errors(
                    "Encountered errors setting up the NPT integrator")

ELSE:
    cdef class VelocityVerletIsotropicNPT(Integrator):
        def __init__(self, *args, **kwargs):
            raise Exception("NPT not compiled in.")


cdef class BrownianDynamics(Integrator):
    """
    Brownian Dynamics integrator.

    """

    def default_params(self):
        return {}

    def valid_keys(self):
        """All parameters that can be set.

        """
        return {}

    def required_keys(self):
        """Parameters that have to be set.

        """
        return {}

    def validate_params(self):
        return True

    def _set_params_in_es_core(self):
        integrate_set_bd()
