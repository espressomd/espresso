# Copyright (C) 2010-2019 The ESPResSo project
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
from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
import numpy as np


@script_interface_register
class MeanVarianceCalculator(ScriptInterfaceHelper):

    """
    Accumulates results from observables.

    Parameters
    ----------
    obs : :class:`espressomd.observables.Observable`
    delta_N : :obj:`int`
        Number of timesteps between subsequent samples for the auto update mechanism.

    Methods
    -------
    update()
        Update the accumulator (get the current values from the observable).

    """
    _so_name = "Accumulators::MeanVarianceCalculator"
    _so_bind_methods = (
        "update",
        "shape",
    )
    _so_creation_policy = "LOCAL"

    def get_mean(self):
        """
        Returns the samples mean values of the respective observable with
        which the accumulator was initialized.
        """
        return np.array(self.call_method("get_mean")).reshape(self.shape())

    def get_variance(self):
        """
        Returns the samples variance for the observable.
        """
        return np.array(self.call_method("get_variance")).reshape(self.shape())


@script_interface_register
class TimeSeries(ScriptInterfaceHelper):

    """
    Records results from observables.

    Parameters
    ----------
    obs : :class:`espressomd.observables.Observable`
    delta_N : :obj:`int`
        Number of timesteps between subsequent samples for the auto update mechanism.

    Methods
    -------
    update()
        Update the accumulator (get the current values from the observable).
    clear()
        Clear the data

    """
    _so_name = "Accumulators::TimeSeries"
    _so_bind_methods = (
        "update",
        "shape",
        "clear"
    )
    _so_creation_policy = "LOCAL"

    def time_series(self):
        """
        Returns the recorded values of the observable.
        """
        return np.array(self.call_method("time_series")).reshape(self.shape())


@script_interface_register
class Correlator(ScriptInterfaceHelper):

    """
    Calculates correlations based on results from observables.

    Parameters
    ----------
    obs1 : :class:`espressomd.observables.Observable`
        The observable :math:`A` to be correlated with :math:`B` (``obs2``).
        If ``obs2`` is omitted, autocorrelation of ``obs1`` is calculated by
        default.

    obs2 : :class:`espressomd.observables.Observable`, optional
        The observable :math:`B` to be correlated with :math:`A` (``obs1``).

    corr_operation : :obj:`str`
        The operation that is performed on :math:`A(t)` and
        :math:`B(t+\\tau)` to obtain :math:`C(\\tau)`. The
        following operations are currently available:

        * ``"scalar_product"``: Scalar product of :math:`A` and
          :math:`B`, i.e., :math:`C=\\sum\\limits_{i} A_i B_i`

        * ``"componentwise_product"``: Componentwise product of
          :math:`A` and :math:`B`, i.e., :math:`C_i = A_i B_i`

        * ``"square_distance_componentwise"``: Each component of
          the correlation vector is the square of the difference
          between the corresponding components of the observables, i.e.,
          :math:`C_i = (A_i-B_i)^2`. Example: when :math:`A` is
          :class:`espressomd.observables.ParticlePositions`, it produces the
          mean square displacement (for each component separately).

        * ``"tensor_product"``: Tensor product of :math:`A` and
          :math:`B`, i.e., :math:`C_{i \\cdot l_B + j} = A_i B_j`
          with :math:`l_B` the length of :math:`B`.

        * ``"complex_conjugate_product"``: assuming that the observables
          consist of a complex and real part
          :math:`A=(A_x+iA_y)`, and :math:`B=(B_x+iB_y)`, this
          operation computes the result :math:`C=(C_x+iC_y)`, as:

          .. math::

                C_x = A_xB_x + A_yB_y\\\\
                C_y = A_yB_x - A_xB_y

        * ``"fcs_acf"``: Fluorescence Correlation Spectroscopy (FCS)
          autocorrelation function, i.e.,

          .. math::

               G_i(\\tau) =
               \\frac{1}{N} \\left< \\exp \\left(
               - \\frac{\\Delta x_i^2(\\tau)}{w_x^2}
               - \\frac{\\Delta y_i^2(\\tau)}{w_y^2}
               - \\frac{\\Delta z_i^2(\\tau)}{w_z^2}
               \\right) \\right>

          where

          .. math::

              \\Delta x_i^2(\\tau) = \\left( x_i(0) - x_i(\\tau) \\right)^2

          is the square displacement of particle
          :math:`i` in the :math:`x` direction, and :math:`w_x`
          is the beam waist of the intensity profile of the
          exciting laser beam,

          .. math::

              W(x,y,z) = I_0 \\exp
              \\left( - \\frac{2x^2}{w_x^2} - \\frac{2y^2}{w_y^2} -
              \\frac{2z^2}{w_z^2} \\right).

          The values of :math:`w_x`, :math:`w_y`, and :math:`w_z`
          are passed to the correlator as ``args``

          The above equations are a
          generalization of the formula presented by Hoefling
          et. al. :cite:`hofling11a`. For more information, see
          references therein. Per each 3 dimensions of the
          observable, one dimension of the correlation output
          is produced. If ``"fcs_acf"`` is used with other observables than
          :class:`espressomd.observables.ParticlePositions`, the physical
          meaning of the result is unclear.

    delta_N : :obj:`int`
        Number of timesteps between subsequent samples for the auto update mechanism.

    tau_max : :obj:`float`
        This is the maximum value of :math:`\\tau` for which the
        correlation should be computed.  Warning: Unless you are using
        the multiple tau correlator, choosing ``tau_max`` of more than
        ``100 * dt`` will result in a huge computational overhead.  In a
        multiple tau correlator with reasonable parameters, ``tau_max``
        can span the entire simulation without too much additional cpu time.

    tau_lin : :obj:`int`
        The number of data-points for which the results are linearly spaced
        in ``tau``. This is a parameter of the multiple tau correlator. If you
        want to use it, make sure that you know how it works. By default, it
        is set equal to ``tau_max`` which results in the trivial linear
        correlator. By setting ``tau_lin < tau_max`` the multiple
        tau correlator is switched on. In many cases, ``tau_lin=16`` is a
        good choice but this may strongly depend on the observables you are
        correlating. For more information, we recommend to read
        ref. :cite:`ramirez10a` or to perform your own tests.

    compress1 : :obj:`str`
        These functions are used to compress the data when
        going to the next level of the multiple tau
        correlator. This is done by producing one value out of two.
        The following compression functions are available:

        * ``"discard2"``: (default value) discard the second value from the time series, use the first value as the result

        * ``"discard1"``: discard the first value from the time series, use the second value as the result

        * ``"linear"``: make a linear combination (average) of the two values

        If only ``compress1`` is specified, then
        the same compression function is used for both
        observables. If both ``compress1`` and ``compress2`` are specified,
        then ``compress1`` is used for ``obs1`` and ``compress2`` for ``obs2``.

        Both ``discard1`` and ``discard2`` are safe for all
        observables but produce poor statistics in the
        tail. For some observables, ``"linear"`` compression
        can be used which makes an average of two
        neighboring values but produces systematic
        errors.  Depending on the observable, the
        systematic error using the ``"linear"`` compression
        can be anything between harmless and disastrous.
        For more information, we recommend to read ref.
        :cite:`ramirez10a` or to perform your own tests.

    compress2 : :obj:`str`, optional
        See ``compress1``.

    args: :obj:`float` of length 3
        Three floats which are passed as arguments to the correlation
        function.  Currently it is only used by ``"fcs_acf"``.
        Other correlation operations will ignore these values.
    """

    _so_name = "Accumulators::Correlator"
    _so_bind_methods = (
        "update",
        "finalize")
    _so_creation_policy = "LOCAL"

    def result(self):
        """
        Returns
        -------

        numpy.ndarray
            The result of the correlation function as a 2d-array.
            The first column contains the values of the lag time tau.
            The second column contains the number of values used to
            perform the averaging of the correlation. Further columns contain
            the values of the correlation function. The number of these columns
            is the dimension of the output of the correlation operation.
        """
        res = np.array(self.call_method("get_correlation"))
        return res.reshape((self.n_result, 2 + self.dim_corr))


@script_interface_register
class AutoUpdateAccumulators(ScriptObjectRegistry):

    """
    Class for handling the auto-update of accumulators used by
    :class:`espressomd.system.System`.

    """
    _so_name = "Accumulators::AutoUpdateAccumulators"
    _so_creation_policy = "LOCAL"

    def add(self, accumulator):
        """
        Adds an accumulator instance to the auto-update list.

        """
        self.call_method("add", object=accumulator)

    def remove(self, accumulator):
        """
        Removes an accumulator from the auto-update list.

        """
        self.call_method("remove", object=accumulator)

    def clear(self):
        """
        Removes all accumulators from the auto-update list.
        """
        self.call_method("clear")
