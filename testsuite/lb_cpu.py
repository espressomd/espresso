from __future__ import print_function
import unittest as ut
import lb_common
import espressomd
import espressomd.lb as lb



if espressomd.has_features("LB"): 
    lb_common.TestLB.lb_class=lb.LBFluid
    lb_common.TestLB.params.update({"mom_prec":1E-9,"mass_prec_per_node":5E-8})
@ut.skipIf(not espressomd.has_features(["LB"]),
           "Features not available, skipping test!")
class TestLBCPU(lb_common.TestLB):
    pass
if __name__ == "__main__":
  ut.main()
