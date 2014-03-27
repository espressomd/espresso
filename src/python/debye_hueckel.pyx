include "myconfig.pxi"
IF ELECTROSTATICS==1:
  def setParams(kappa, rCut):
      if rCut<0:
          raise ValueError("rCut must be > 0")
      dh_set_params(kappa, rCut)
  
  def setRcut(rCut):
      if rCut<0:
          raise ValueError("rCut must be > 0")
      dh_set_params(dh_params.kappa, rCut)
  
  def setKappa(kappa):
      dh_set_params(kappa, dh_params.r_cut)
  
  def getRcut():
      return dh_params.r_cut
  
  def getKappa():
      return dh_params.kappa
  
