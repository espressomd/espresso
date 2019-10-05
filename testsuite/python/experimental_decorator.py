#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd
from espressomd.utils import requires_experimental_features

import unittest as ut


@requires_experimental_features("because")
class A:

    def __init__(self, *args, **kwargs):
        pass


class RequireExperimental(ut.TestCase):

    def test(self):

        if espressomd.has_features("EXPERIMENTAL_FEATURES"):
            x = A(1, 3.0)
        else:
            with self.assertRaisesRegexp(Exception, "because"):
                x = A(1, 3.0)


if __name__ == "__main__":
    ut.main()
