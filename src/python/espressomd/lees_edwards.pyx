from .script_interface import ScriptInterfaceHelper, script_interface_register
from .utils import handle_errors, to_str

from globals cimport *

lees_edwards_type_dict = {'off': LEES_EDWARDS_PROTOCOL_OFF,
                          'step': LEES_EDWARDS_PROTOCOL_STEP,
                          'steady_shear': LEES_EDWARDS_PROTOCOL_STEADY_SHEAR,
                          'oscillatory_shear': LEES_EDWARDS_PROTOCOL_OSC_SHEAR}


@script_interface_register
class LeesEdwards(ScriptInterfaceHelper):

    """Inteface to the Lees-Edwards binding.

       See documentation.

    """

    _so_name = "LeesEdwards::LeesEdwards"

    def __init__(self, *args, **kwargs):
        """
        This class should not be instanced by the user. Instead, use
        the :attr:`espressomd.system.System.lees_edwards` attribute
        of the system class to access the Lees-Edwards binding.

        Use :meth:`espressomd.lees_edwards.LeesEdwards.set_params`
        to change the parameters of the Lees-Edwards boundary.

        """
        # If no mode is specified at construction, use off.
        if "type" not in kwargs:
            kwargs["type"] = "off"
        super(type(self), self).__init__()
        self.set_params(**kwargs)

    def validate(self):
        """Validates the parameters of the Lees-Edwards boundary.

           This is called automatically on parameter change

        """
        return self.call_method("validate")

    # Do not allow setting of individual attributes
    def __setattr__(self, *args, **kwargs):
        raise Exception(
            "Please set all parameters at once via lees_edwards.set_params()")

    # Override to call validate after parameter update
    def set_params(self, **kwargs):
        """
        Set the parameters for the Lees-Edwards boundary.

        See the documentation.


        Parameters
        ----------
        type : One of "off", "step", "steady_shear", "oscillatory_shear".
               Lees-Edwards protocol               

        offset : :obj:`float`
               Fixed offset for "step" protocol

        velocity : :obj:`float`
               Fixed velocity for "steady_shear" protocol

        frequency : :obj:`float`
               Frequency for "oscillatory_shear" protocol

        amplitude : :obj:`float`
               Amplitude for "oscillatory_shear" protocol

        sheardir : :obj:`int`
               Direction of shear flow

        shearplanenormal : :obj:`int`
               Direction normal to the shear flow

        """
        kwargs["type"] = lees_edwards_type_dict[kwargs["type"]]

        if not ("type" in kwargs):
            raise Exception(
                "Lees-Edwards protocol type must be specified via the type keyword argument")

        if kwargs['type'] == 'step':
            if 'offset' not in kwargs:
                raise Exception('No offset given for step strain')
        elif kwargs['type'] == 'steady_shear':
            if 'velocity' not in kwargs:
                raise Exception('No velocity given for steady shear')
        elif kwargs['type'] == 'oscillatory_shear':
            if 'amplitude' not in kwargs or 'frequency' not in kwargs:
                raise Exception('No amplitude or frequency given for oscillatory shear')
        else:
            raise Exception('Lees-Edwards protocol '+ str(kwargs['type']) +' unknown')

        if 'sheardir' not in kwargs:
            kwargs['sheardir'] = 0

        if 'shearplanenormal' not in kwargs:
            kwargs['shearplanenormal'] = 1

        kwargs['time0'] = sim_time

        super(type(self), self).set_params(**kwargs)

        self.validate()

        handle_errors("Validation of Lees-Edwards parameters failed")

    def get_parameter(self, name):
        """Gets a single parameter from the Lees-Edwards boundary."""

        res = super(type(self), self).get_parameter(name)
        return self._convert_param(name, res)

    def get_params(self):
        """Returns the parameters of the Lees-Edwards boundary as a dict."""

        res = super(type(self), self).get_params()
        for k in res.keys():
            res[k] = self._convert_param(k, res[k])

        # Filter key-value pairs according to active type
        return {k: res[k] for k in self._params_for_mode(res["type"])}

    def _convert_param(self, name, value):
        """Handles type conversion core -> python"""
        # Py3: Cast from binary to normal string. Don't understand, why a
        # binary string can even occur, here, but it does.
        name = to_str(name)
        # Convert int mode parameter to string
        res = value
        if name == "type":
            res = self._str_type(value)

        return res

    def _params_for_mode(self, type):
        """The parameter names expected for a given collision mode

        """
        if type == "off":
            return ("type", "sheardir", "shearplanenormal")
        if type == "step":
            return ("type", "offset", "sheardir", "shearplanenormal")
        if type == "steady_shear":
            return ("type", "velocity", "sheardir", "shearplanenormal")
        if type == "oscillatory":
            return ("type", "amplitude", "frequency", "sheardir", "shearplanenormal")

        raise Exception("Type not handled: " + type.__str__())

    def _str_type(self, int_type):
        """String type name from int ones provided by the core

        """
        for key in lees_edwards_type_dict:
            if lees_edwards_type_dict[key] == int_type:
                return key
        raise Exception("Unknown integer Lees-Edwards protocol type %d" % int_type)
