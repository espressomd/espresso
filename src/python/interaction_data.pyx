# Non-bonded interactions

cdef class NonBondedInteractionHandle:
  """Provides access to Non-bonded interaction parameters for a specific pair of particle types """
  cdef int type1
  cdef int type2
  cdef IA_parameters* params

  def __init__(self, _type1, _type2):
    """Takes two particle types as argument"""
    self.type1=_type1
    self.type2=_type2
    self.update()
  
  def update(self):
    """Fetches the current interaction parameters from the Espresso core"""
    self.params =  get_ia_param(self.type1,self.type2) 

  property lennardJones:
    """Lennard Jones parameters"""
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




cdef class NonBondedInteractions:
  """Access to non-bonded interaction parameters via [i,j], where i,j are particle 
  types. Returns NonBondedInteractionHandle."""
  def __getitem__(self,key):
    if not isinstance(key,tuple):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    if len(key) != 2 or (not isinstance(key[0],int)) or (not isinstance(key[1],int)):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    return NonBondedInteractionHandle(key[0],key[1])





# Bonded interactions

cdef class BondedInteraction(object):
  def __init__(self, *args, **kwargs):
    """Either called with an interaction id, in which case, the interaction will represent
       the bonded interaction as it is defined in Espresso core
       Or called with keyword arguments describing a new interaction."""
    # Interaction id as argument
    if len(args)==1 and isinstance(args[0],int):
      bondId=args[0]
      # Check, if the bond in Espresso core is really defined as a FENE bond
      if bonded_ia_params[bondId].type != self.typeNumber():
        raise Exception("The bond with this id is not defined as a "+self.typeName()+" bond in the Espresso core.")
      
      # Load the parameters currently set in the Espresso core
      self._params=self._getParamsFromEsCore()
      self._bondId=bondId
    
    # Or have we been called with keyword args describing the interaction
    elif len(args)==0:
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__)
      
      # Validation of parameters
      self.validateParams()
      
    else: 
      raise Exception("The constructor has to be called either with a bond id (as interger), or with a set of keyword arguments describing a new interaction")

   
      
   


  def isValid(self):
    """Check, if the data stored in the instance still matches what is in Espresso"""
    # Check if the bond type in Espresso still matches the bond type saved in this class
    if bonded_ia_params[self._bondId].type != self.typeNumber():
      return False

    # check, if the bond parameters saved in the class still match those saved in Espresso
    tempParams =self._getParamsFromEsCore()
    if self._params != tempParams:
      return False
    
    # If we're still here, the instance is valid
    return True
 
 
  property params:
    def __get__(self):
      return self._params

    def __set__(self,p):
      # Check, if any key was passed, which is not known
      for k in p.keys:
        if k not in self.validKeys():
          raise ValueError("Only the following keys are supported: "+self.validKeys().__str__)
      
      # Initialize default values
      self.setDefaultParams()
      # Put in values given by the user
      self._params.update(p)

  def validateParameters(self):
    return True

  def _getParamsFromEsCore(self):
    raise Exception("Subclasses of BondedInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of BondedInteraction must define the setParamsFromEsCore() method.")
  
  def setDefaultParams(self):
    raise Exception("Subclasses of BondedInteraction must define the setDefaultParams() method.")

  def typeNumber(self): 
    raise Exception("Subclasses of BondedInteraction must define the typeNumber() method.")
   
  
  def typeName(self): 
    raise Exception("Subclasses of BondedInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of BondedInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of BondedInteraction must define the requiredKeys() method.")





# Fene bond

class FeneBond(BondedInteraction):

  def typeNumber(self):
    return 0

  def typeName(self): 
    return "FENE"

  def validKeys(self):
    return "k","d_r_max","r_0"

  def requiredKeys(self): 
    return "k","d_r_max"

  def setDefaultParams(self):
    self._params = {"r_0":0.} 
    # Everything else has to be supplied by the user, anyway

  def _getParamsFromEsCore(self):
    return \
      {"k":bonded_ia_params[self._bondId].p.fene.k,\
       "d_r_max":bonded_ia_params[self._bondId].p.fene.drmax,\
       "r_0":bonded_ia_params[self._bondId].p.fene.r0}

  def _setParamsInEsCore(self):
   fene_set_params(self._bondId,self._params["k"],self._params["d_r_max"],self.params["r_0"])

class HarmonicBond(BondedInteraction):
  pass
    


bondedInteractionClasses = {0:FeneBond, 1:HarmonicBond}






class BondedInteractions:
  """Represents the non-bonded interactions. Individual interactions can be accessed using
  NonBondedInteractions[i], where i is the bond id. Will return an instance o
  BondedInteractionHandle"""
  def __getItem__(self, key):
    if not isinstance(key,int):
      raise ValueError("Index to BondedInteractions[] hast to ba an integer referring to a bond id")

    # Find out the type of the interaction from Espresso
    bondType = bonded_ia_params[key].type

    # Check if the bonded interaction exists in Espresso core
    if bondType == -1:
      raise ValueError("The bonded interaction with the id "+str(key)+" is not yet defined.")

    # Find the appropriate class representing such a bond
    bondClass =bondedInteractionClasses[key]

    # And return an instance of it, which refers to the bonded interaction id in Espresso
    return bondClass(key)
     
  def __setitem__(self,key,value):
    # Validate arguments
   
    # type of key must be int
    if not isinstance(key,int):
      raise ValueError("Index to BondedInteractions[] has to ba an integer referring to a bond id")

    # Value must be subclass off BondedInteraction
    if not isinstance(value,BondedInteraction):
      raise ValueError("Only subclasses of BondedInteraction can be assigned.")

    # Save the bond id in the BondedInteraction instance
    value._bondId=key

    # Set the parameters of the BondedInteraction instance in the Es core
    value._setParamsInEsCore()






