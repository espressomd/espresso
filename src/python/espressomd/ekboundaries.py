from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper, script_interface_register
import espressomd.lbboundaries

if any(i in espressomd.code_info.features() for i in ["LB_BOUNDARIES", "LB_BOUNDARIES_GPU"]):
    @script_interface_register
    class EKBoundaries(espressomd.lbboundaries.LBBoundaries):

        """
        Creates a set of electrokinetics boundaries.

        """
        pass

    @script_interface_register
    class EKBoundary(espressomd.lbboundaries.LBBoundary):

        """
        Creates a EK boundary.

        """
        pass
