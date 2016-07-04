#
# Copyright (C) 2014,2015,2016 The ESPResSo project
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
cdef extern from "tcl.h":
    cdef struct Tcl_Interp:
        pass

    Tcl_Interp * Tcl_CreateInterp()
    void Tcl_DeleteInterp(Tcl_Interp *)
    int Tcl_Eval(Tcl_Interp * interp, char * script)
    char * Tcl_GetStringResult(Tcl_Interp * interp)

cdef extern from "tcl/initialize_interpreter.hpp":
    int tcl_appinit(Tcl_Interp * interp)

# Define Tcl interpreter
cdef class TclInterpreter:
    cdef Tcl_Interp * interp

    def __init__(self):
        self.interp = Tcl_CreateInterp()
        self.eval('global argv; set argv ""')
        self.eval('set tcl_interactive 0')
        tcl_appinit(self.interp)

    def eval(self, string):
        cdef const char * tclresult
        result = Tcl_Eval(self.interp, string)
        tclresult = Tcl_GetStringResult(self.interp)
        if result:
            raise Exception("Tcl reports an error", tclresult)

        return tclresult

    def __del__(self):
        Tcl_DeleteInterp(self.interp)
