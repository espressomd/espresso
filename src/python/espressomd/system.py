from __future__ import print_function
from espressomd.highlander import highlander
from espressomd import _system


@highlander
class System(_system.System):

    def __init__(self):
        _system.System.__init__(self)

    def __setattr__(self, name, value):
        if hasattr(self, name):
            super(System, self).__setattr__(name, value)
        else:
            raise AttributeError(
                "System does not have the attribute " + name + "."
                + "\nIf you know what you're doing, use "
                + "system.create_attr()"
            )

    def create_attr(self, *args, **kwargs):
        """
        Circumvents the __setattr__ lock. Allows to create and set new
        attributes to System instances.
        For *args, it initializes an attribute with None value, if the
        attribute does not already exist.
        For **kwargs, it simply calls super().__setattr__.
        """
        for arg in args:
            try: name = str(arg)
            except ValueError:
                print("Please pass either **kwargs or string *args to"
                      + "create_attr()"
                      )
                continue

            if hasattr(self, name):
                print("Attribute " + name + " already exists.")
            else:
                super(System, self).__setattr__(name, None)

        for name, value in list(kwargs.items()):
            super(System, self).__setattr__(name, value)
