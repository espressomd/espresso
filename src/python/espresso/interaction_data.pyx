cdef class NonBondedInteractionHandle:
  cdef int type1
  cdef int type2
  cdef IA_parameters* params

  def __init__(self, _type1, _type2):
    self.type1=_type1
    self.type2=_type2
    self.update()
  
  def update(self):
    self.params =  get_ia_param(self.type1,self.type2) 

  property lennardJones:
    def __set__(self, value):
      self.update()

      if "eps" in value:
        self.params[0].LJ_eps =value["eps"]
        del value["eps"]
      
      if "sigma" in value:
        self.params[0].LJ_sig =value["sigma"]
        del value["sigma"]

      if "cut" in value:
        self.params[0].LJ_cut =value["cut"]
        del value["cut"]
      
      if "shift" in value:
        self.params[0].LJ_shift =value["shift"]
        del value["shift"]

      if "offset" in value:
        self.params[0].LJ_offset =value["offset"]
        del value["offset"]

      if "ljcap" in value:
        ljcap =value["ljcap"]
        del value["ljcap"]

# Currently not working, since overwritten by c:
#      if "capradius" in value:
#        self.params[0].LJ_capradius =value["capradius"]
#        del value["capradius"]

      if "min" in value:
        self.params[0].LJ_min =value["min"]
        del value["min"]

      if (len(value) >0):
        raise Exception("Unsupported parameters: " +value.__str__())

      lennard_jones_set_params(self.type1, self.type2, 
         self.params[0].LJ_eps,self.params[0].LJ_sig, self.params[0].LJ_cut,
         self.params[0].LJ_shift,self.params[0].LJ_offset, self.params[0].LJ_capradius,
         self.params[0].LJ_min)

      ljforcecap_set_params(ljcap)
      
      self.update()


    def __get__(self):
      self.update()
      return {"eps":self.params[0].LJ_eps,
              "sigma":self.params[0].LJ_sig,
              "cut":self.params[0].LJ_cut,
              "shift":self.params[0].LJ_shift,
              "offset":self.params[0].LJ_offset,
              "capradius":self.params[0].LJ_capradius,
              "min":self.params[0].LJ_min}


#class InteractionList:
#  def __getItem__(key):
#    return NonBondedInteractionHandle(0,0)


class InteractionList:
  def __getitem__(self,key):
    return NonBondedInteractionHandle(key[0],key[1])

