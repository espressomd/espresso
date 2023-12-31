#
# Copyright (C) 2023 The ESPResSo project
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


def update_file(filepath, function):
    with open(filepath, "r+") as f:
        content = function(f.read())
        f.seek(0)
        f.truncate()
        f.write(content)


def update_espressomd(content):
    token = ".. automodule:: espressomd.propagation"
    assert token in content
    assert content.count(token) == 1
    content = content.replace(token, token + "\n   :member-order: bysource", 1)
    return content


update_file("espressomd.rst", update_espressomd)
