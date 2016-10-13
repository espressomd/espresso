from ...script_interface import PScriptInterface
import sys


class H5md():
    def __init__(self, **kwargs):
        self.valid_params = ['filename',]
        if not ('filename' in kwargs):
            raise ValueError("'filename' parameter missing.")
        self.what = {'write_pos': 1<<0,
                     'write_vel': 1<<1,
                     'write_force': 1<<2,
                     'write_type': 1<<3,
                     'write_mass': 1<<4}
        self.valid_params.append(self.what.keys())
        self.what_bin = 0
        for i,j in kwargs.iteritems():
            if (i in self.what.keys()):
                if (isinstance(j, bool)):
                    if j:
                        self.what_bin += self.what[i]
                else:
                    raise ValueError("{} has to be a bool value.".format(i))
            elif (i not in self.valid_params):
                raise ValueError("Unknown parameter {} for H5MD writer.".format(i))
        self.h5md_instance = PScriptInterface("ScriptInterface::Writer::H5mdScript")
        self.h5md_instance.set_params(filename=kwargs['filename'], what=self.what_bin, 
                                      scriptname=sys.argv[0])
        self.h5md_instance.call_method("init_file")

    def _print_params(self):
        return self.h5md_instance.get_params()


    def write(self):
        self.h5md_instance.call_method("write")
    

    def close(self):
        self.h5md_instance.call_method("close")
