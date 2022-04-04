# Copyright (C) 2021 The ESPResSo project
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

import espressomd
import unittest as ut
import unittest_decorators as utx
import numpy as np


class Tests(ut.TestCase):

    def get_skip_reason(self, decorator):
        try:
            decorator(lambda: None)()
        except ut.case.SkipTest as err:
            return err.args

    def test_version_requirements(self):
        fun = utx.skipIfUnmetModuleVersionRequirement
        decorator = fun('numpy___', '>' + np.__version__)
        args = self.get_skip_reason(decorator)
        self.assertEqual(
            args, ('Skipping test: missing python module numpy___',))
        decorator = fun('numpy', '>1' + np.__version__)
        args = self.get_skip_reason(decorator)
        self.assertEqual(
            args, ('Skipping test: version requirement not met for module numpy',))
        decorator = fun('numpy', '==' + np.__version__)
        args = self.get_skip_reason(decorator)
        self.assertIsNone(args)

    def test_missing_modules(self):
        err_msg = 'Skipping test: missing python '
        decorator = utx.skipIfMissingModules('numpy___')
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, (err_msg + 'module numpy___',))
        decorator = utx.skipIfMissingModules('numpy___', 'scipy___')
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, (err_msg + 'modules numpy___, scipy___',))
        decorator = utx.skipIfMissingModules(['numpy___', 'scipy___'])
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, (err_msg + 'modules numpy___, scipy___',))
        decorator = utx.skipIfMissingModules(['espressomd'])
        args = self.get_skip_reason(decorator)
        self.assertIsNone(args)

    def test_missing_gpu(self):
        espressomd.gpu_available = lambda: False
        decorator = utx.skipIfMissingGPU()
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, ('Skipping test: no GPU available',))
        espressomd.gpu_available = lambda: True
        decorator = utx.skipIfMissingGPU()
        args = self.get_skip_reason(decorator)
        self.assertIsNone(args)

    def test_missing_features(self):
        err_msg = 'Skipping test: missing features '
        espressomd.has_features = lambda x: False
        decorator = utx.skipIfMissingFeatures('UNKNOWN')
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, (err_msg + 'UNKNOWN',))
        decorator = utx.skipIfMissingFeatures(['UNKNOWN1', 'UNKNOWN2'])
        args = self.get_skip_reason(decorator)
        self.assertEqual(args, (err_msg + 'UNKNOWN1, UNKNOWN2',))
        espressomd.has_features = lambda x: True
        decorator = utx.skipIfMissingFeatures('KNOWN')
        args = self.get_skip_reason(decorator)
        self.assertIsNone(args)


if __name__ == '__main__':
    ut.main()
