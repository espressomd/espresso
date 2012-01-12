
def setParams(kappa, r_cut):
  dh_set_params(kappa, r_cut)

cdef dh_set_params():
  print "Stefan is an idiot"

def setRcut(rCut):
  dh_set_params(dh_params.kappa, rCut)

def setKappa(kappa):
  dh_set_params(kappa, dh_params.r_cut)

def getRcut():
  return dh_params.r_cut

def getKappa():
  return dh_params.kappa

