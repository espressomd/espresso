#
# Copyright (C) 2013-2022 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .code_features import assert_features
import signal


@script_interface_register
class IntegratorHandle(ScriptInterfaceHelper):
    """
    Provide access to the currently active integrator.

    """
    _so_name = "Integrators::IntegratorHandle"
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        if "sip" not in kwargs and "integrator" not in kwargs:
            kwargs["integrator"] = VelocityVerlet()
        super().__init__(**kwargs)

    def __str__(self):
        return f'{self.__class__.__name__}({self.integrator.__class__.__name__})'

    def run(self, *args, **kwargs):
        """
        Run the integrator.

        """
        return self.integrator.run(*args, **kwargs)

    def set_steepest_descent(self, **kwargs):
        """
        Set the integration method to steepest descent
        (:class:`SteepestDescent`).

        """
        self.integrator = SteepestDescent(**kwargs)

    def set_vv(self):
        """
        Set the integration method to velocity Verlet, which is suitable for
        simulations in the NVT ensemble (:class:`VelocityVerlet`).

        """
        self.integrator = VelocityVerlet()

    def set_nvt(self):
        """
        Set the integration method to velocity Verlet, which is suitable for
        simulations in the NVT ensemble (:class:`VelocityVerlet`).

        """
        self.integrator = VelocityVerlet()

    def set_isotropic_npt(self, **kwargs):
        """
        Set the integration method to a modified velocity Verlet designed for
        simulations in the NpT ensemble (:class:`VelocityVerletIsotropicNPT`).

        """
        self.integrator = VelocityVerletIsotropicNPT(**kwargs)

    def set_brownian_dynamics(self):
        """
        Set the integration method to BD.

        """
        self.integrator = BrownianDynamics()

    def set_stokesian_dynamics(self, **kwargs):
        """
        Set the integration method to Stokesian Dynamics (:class:`StokesianDynamics`).

        """
        self.integrator = StokesianDynamics(**kwargs)


class Integrator(ScriptInterfaceHelper):
    """
    Integrator class.

    """

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
        status = self.call_method(
            "integrate",
            with_nogil=True,
            steps=steps,
            recalc_forces=recalc_forces,
            reuse_forces=reuse_forces)
        self.handle_sigint(status)

    def handle_sigint(self, error_code):
        if error_code == -2:
            signal.raise_signal(signal.Signals.SIGINT)


@script_interface_register
class SteepestDescent(Integrator):
    """
    Steepest descent algorithm for energy minimization.

    Particles located at :math:`\\vec{r}_i` at integration step :math:`i` and
    experiencing a potential :math:`\\mathcal{H}(\\vec{r}_i)` are displaced
    according to the equation:

    :math:`\\vec{r}_{i+1} = \\vec{r}_i - \\gamma\\nabla\\mathcal{H}(\\vec{r}_i)`

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

    """
    _so_name = "Integrators::SteepestDescent"
    _so_creation_policy = "GLOBAL"

    def run(self, steps=1):
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
        integrated = self.call_method(
            "integrate", with_nogil=True, steps=steps)
        self.handle_sigint(integrated)

        return integrated


@script_interface_register
class VelocityVerlet(Integrator):
    """
    Velocity Verlet integrator, suitable for simulations in the NVT ensemble.

    """
    _so_name = "Integrators::VelocityVerlet"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class VelocityVerletIsotropicNPT(Integrator):
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
        them to ``True``. Default is all ``True``.
    cubic_box : :obj:`bool`, optional
        If ``True``, a cubic box is assumed and the value of ``direction``
        will be ignored when rescaling the box. This is required e.g. for
        electrostatics and magnetostatics. Default is ``False``.

    """
    _so_name = "Integrators::VelocityVerletIsoNPT"
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        assert_features("NPT")
        super().__init__(**kwargs)


@script_interface_register
class BrownianDynamics(Integrator):
    """
    Brownian Dynamics integrator.

    """
    _so_name = "Integrators::BrownianDynamics"
    _so_creation_policy = "GLOBAL"


@script_interface_register
class StokesianDynamics(Integrator):
    """
    Stokesian Dynamics integrator.

    Parameters
    ----------
    viscosity : :obj:`float`
        Bulk viscosity.
    radii : :obj:`dict`
        Dictionary that maps particle types to radii.
    approximation_method : :obj:`str`, optional, {'ft', 'fts'}
        Chooses the method of the mobility approximation.
        ``'fts'`` is more accurate. Default is ``'fts'``.
    self_mobility : :obj:`bool`, optional
        Switches off or on the mobility terms for single particles. Default
        is ``True``.
    pair_mobility : :obj:`bool`, optional
        Switches off or on the hydrodynamic interactions between particles.
        Default is ``True``.

    """
    _so_name = "Integrators::StokesianDynamics"
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        assert_features("STOKESIAN_DYNAMICS")
        if kwargs.get("lubrication", False):
            raise NotImplementedError(
                "Stokesian Dynamics lubrication is not available yet")
        super().__init__(**kwargs)
