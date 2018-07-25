from __future__ import print_function
import unittest as ut
import lb_common
import espressomd
import espressomd.lb as lb



if espressomd.has_features("LB_GPU"): 
    lb_common.TestLB.lb_class=lb.LBFluidGPU
    lb_common.TestLB.params.update({"mom_prec":1E-3,"mass_prec_per_node":1E-5})
@ut.skipIf(not espressomd.has_features(["LB_GPU"]),
           "Features not available, skipping test!")
@ut.skipIf(espressomd.has_features("SHANCHEN"), "Test not compatible with Shan Chen")
class TestLBGPU(lb_common.TestLB):
    pass
if __name__ == "__main__":
  ut.main()
