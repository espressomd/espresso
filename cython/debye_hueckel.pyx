
def setRcut(rcut):
  dh_set_params(dh_params.kappa, rcut)

def getRcut():
  return dh_params.rcut
