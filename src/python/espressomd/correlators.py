from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
import numpy as np


@script_interface_register
class Correlator(ScriptInterfaceHelper):
    """
    Calculates correlations based on results from observables.

    Parameters
    ----------
    obs1, obs2 : Instances of :class:`espressomd.observables.Observable`.
                 The observables A and B that are to be correlated. If `obs2`
                 is omitted, autocorrelation of `obs1` is calculated by
                 default.
    corr_operation : :obj:`str`
                     The operation that is performed on :math:`A(t)` and
                     :math:`B(t+\\tau)` to obtain :math:`C(\\tau)`. The
                     following operations are currently available:

                         * `scalar_product`: Scalar product of :math:`A` and
                           :math:`B`, i.e., :math:`C=\sum\limits_{i} A_i B_i`

                         * `componentwise_product`: Componentwise product of
                           :math:`A` and :math:`B`, i.e., :math:`C_i = A_i B_i`

                         * `square_distance_componentwise`: Each component of
                           the correlation vector is the square of the difference
                           between the corresponding components of the
                           observables, i.E., :math:`C_i = (A_i-B_i)^2`. Example:
                           when :math:`A` is `ParticlePositions`, it produces the
                           mean square displacement (for each component
                           separately).

                         * `tensor_product`: Tensor product of :math:`A` and
                           :math:`B`, i.e., :math:`C_{i \\cdot l_B + j} = A_i B_j`
                           with :math:`l_B` the length of :math:`B`.

                         * `complex_conjugate_product`:

                           .. todo:: `Complex conjugate product must be defined.`

                         * `fcs_acf`:

                           .. todo:: `How to set the extra args in py?`

                           Fluorescence Correlation Spectroscopy (FCS)
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

                           The equations are a
                           generalization of the formula presented by Hoefling
                           et. al. :cite:`hofling11a`. For more information, see
                           references therein. Per each 3 dimensions of the
                           observable, one dimension of the correlation output
                           is produced. If `fcs_acf` is used with other
                           observables than `ParticlePositions`, the physical
                           meaning of the result is unclear.
    dt : :obj:`float`
         The time interval of sampling data points. When autoupdate is used,
         `dt` has to be a multiple of timestep. It is also used to produce time
         axis in real units. Warning: if `dt` is close to the timestep,
         autoupdate is strongly recommended. Otherwise cpu time is wasted on
         passing the control between the script and kernel.

    tau_max : :obj:`float`
              This is the maximum value of :math:`\tau` for which the
              correlation should be computed.  Warning: Unless you are using
              the multiple tau correlator, choosing `tau_max` of more than
              100`dt` will result in a huge computational overhead.  In a
              multiple tau correlator with reasonable parameters, `tau_max`
              can span the entire simulation without too much additional cpu
              time.

    tau_lin : :obj:`int`
              The number of data-points for which the results are linearly spaced
              in `tau`. This is a parameter of the multiple tau correlator. If you
              want to use it, make sure that you know how it works. By default, it
              is set equal to `tau_max` which results in the trivial linear
              correlator. By setting `tau_lin` < `tau_max` the multiple
              tau correlator is switched on. In many cases, `tau_lin`=16 is a
              good choice but this may strongly depend on the observables you are
              correlating. For more information, we recommend to read
              Ref. :cite:`ramirez10a` or to perform your own tests.

    compress1 and compress2 : :obj:`str`
                              Are functions used to compress the data when
                              going to the next level of the multiple tau
                              correlator. Different compression functions for
                              different observables can be specified if
                              desired, otherwise the same function is used for
                              both.  Default is `discard` which takes one of
                              the observable values and discards the other one.
                              This is safe for all observables but produces
                              poor statistics in the tail. For some
                              observables, linear compression can be used
                              which makes an average of two neighboring values
                              but produces systematic errors.  Depending on the
                              observable, the systematic error can be anything
                              between harmless and disastrous. For more
                              information, we recommend to read
                              Ref. :cite:`ramirez10a` or to perform your own
                              tests.

    """

    _so_name = "Correlators::Correlator"
    _so_bind_methods = (
        "update",
        "auto_update",
        "finalize",
        "dim_corr",
        "n_results",
        "hierarchy_depth")
    _so_creation_policy = "LOCAL"

    def result(self):
        res = np.array(self.call_method("get_correlation"))
        res = res.reshape((self.n_results(), 2 + self.dim_corr()))
        return res


@script_interface_register
class AutoUpdateCorrelators(ScriptInterfaceHelper):
    _so_name = "Correlators::AutoUpdateCorrelators"
    _so_creation_policy = "LOCAL"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Correlator):
                correlator = args[0]
            else:
                raise TypeError(
                    "Either a Correlator object or key-value pairs for the parameters of a Correlator object need to be passed.")
        else:
            correlator = Correlator(**kwargs)
        self.call_method("add", object=correlator)
        return correlator

    def remove(self, Correlator):
        self.call_method("remove", object=Correlator)
