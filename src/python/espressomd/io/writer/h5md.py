# -*- coding: utf-8 -*-
"""Interface module for the H5md core implementation.
"""


import sys
from six import iteritems
from ...script_interface import PScriptInterface # pylint: disable=import-error


class H5md(object):
    """H5md file object for accessing to the H5MD core
       implementation via the PScriptInterface.
    """

    def __init__(self, **kwargs):
        """Initialize a H5md object.

        Parameters
        ----------
        filename : str
                   Name of the trajectory file.
        write_pos : bool, optional
                    If positions should be written.
        write_vel : bool, optional
                    If velocities should be written.
        write_force : bool, optional
                      If forces should be written.
        write_type : bool, optional
                     If types should be written.
        write_mass : bool, optional
                     If masses should be written.
        write_charge : bool, opional
                       If charges should be written.
        write_ordered: bool, optional
                       If particle properties should be ordered according to ids.
        """
        self.valid_params = ['filename', "write_ordered"]
        if 'filename' not in kwargs:
            raise ValueError("'filename' parameter missing.")
        self.what = {'write_pos': 1<<0,
                     'write_vel': 1<<1,
                     'write_force': 1<<2,
                     'write_type': 1<<3,
                     'write_mass': 1<<4,
                     'write_charge': 1<<5}
        self.valid_params.append(self.what.keys())
        self.what_bin = 0
        for i, j in iteritems(kwargs):
            if i in self.what.keys():
                if isinstance(j, bool):
                    if j:
                        self.what_bin += self.what[i]
                else:
                    raise ValueError("{} has to be a bool value.".format(i))
            elif i not in self.valid_params:
                raise ValueError("Unknown parameter {} for H5MD writer.".format(i))
        
        write_ordered_default=True
        self.write_ordered=kwargs.get('write_ordered', write_ordered_default)
        
        self.h5md_instance = PScriptInterface("ScriptInterface::Writer::H5mdScript")
        self.h5md_instance.set_params(filename=kwargs['filename'], what=self.what_bin,
                                      scriptname=sys.argv[0], write_ordered=self.write_ordered)
        self.h5md_instance.call_method("init_file")


    def get_params(self):
        """Get the parameters from the scriptinterface.
        """
        return self.h5md_instance.get_params()


    def write(self):
        """Calls the H5md write method.
        """
        self.h5md_instance.call_method("write")

    def flush(self):
        """Calls the H5md flush method.
        """
        self.h5md_instance.call_method("flush")

    def close(self):
        """Close the H5md file.
        """
        self.h5md_instance.call_method("close")
