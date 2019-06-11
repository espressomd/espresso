# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import tempfile
import subprocess
from sys import argv


class Defines(object):

    def __init__(self, compiler, flags=[]):
        self._compiler = compiler
        self._flags = flags + ["-E", "-dM"]

        build_in = self._build_in_defs()
        self._buildin = set(build_in)

    def _run_compiler(self, filename):
        raw_output = subprocess.check_output(
            [self._compiler] + self._flags + [filename], stderr=subprocess.STDOUT)
        return raw_output.decode('ascii').splitlines()

    def _remove_define(self, line):
        return line[len("#define"):].strip()

    def _get_defs(self, filename):
        lines = self._run_compiler(filename)
        return map(self._remove_define, lines)

    def _build_in_defs(self):
        with tempfile.NamedTemporaryFile(delete=True, suffix='.cpp') as empty_file:
            return self._get_defs(empty_file.name)

    def build_in_defines(self):
        return self._buildin

    def defines(self, filename, include_build_in=False):
        all_defs = set(self._get_defs(filename))

        if include_build_in:
            return all_defs
        else:
            return all_defs - self._buildin

if __name__ == "__main__":
    compiler = argv[1]
    filename = argv[2]
    flags = argv[3:]

    parser = Defines(compiler, flags)
    map(print, parser.defines(filename))
