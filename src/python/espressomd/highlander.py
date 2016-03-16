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
from functools import update_wrapper


class ThereCanOnlyBeOne(BaseException):

    def __init__(self, cls):
        self._cls = cls

    def __str__(self):
        return "There can only be one instance of '{}' at any time.".format(self._cls)


def highlander(klass):
    klass.highlander_created = False

    def cls_init(self, *args, **kwargs):
        "__init__ method by the highlander decorator."
        if self.__class__.highlander_created:
            raise ThereCanOnlyBeOne(self.__class__)
        self.__class__.highlander_created = True

    def cls_init_call_orig(self, *args, **kwargs):
        if self.__class__.highlander_created:
            raise ThereCanOnlyBeOne(self.__class__)
        self.__class__.highlander_created = True
        self.__class__.__init_orig__(self, *args, **kwargs)

    # override the __init__ method of the class to store the bool
    # "highlander_created"
    if hasattr(klass, '__init__'):
        klass.__init_orig__ = klass.__init__
        klass.__init__ = cls_init_call_orig
        update_wrapper(cls_init_call_orig, klass.__init_orig__)
    else:
        klass.__init__ = cls_init

    # override the __del__ method of the class
    def cls_del(self):
        "__del__ method by the highlander decorator."
        self.__class__.highlander_created = False

    def cls_del_call_orig(self):
        cls_del(self)
        self.__class__.__del_orig__(self)

    if hasattr(klass, '__del__'):
        klass.__del_orig__ = klass.__del__
        klass.__del__ = cls_del_call_orig
        update_wrapper(cls_del_call_orig, klass.__del_orig__)
    else:
        klass.__del__ = cls_del

    return klass
