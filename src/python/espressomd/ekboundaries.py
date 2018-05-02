from __future__ import print_function, absolute_import

import espressomd
import espressomd.code_info
import espressomd.lbboundaries

if "EK_BOUNDARIES" in espressomd.code_info.features():
    class EKBoundaries(espressomd.lbboundaries.LBBoundaries):
        """
        Creates a set of electrokinetics boundaries.

        """
        pass

    class EKBoundary(espressomd.lbboundaries.LBBoundary):
        """
        Creates a EK boundary.
    
        """
        pass
