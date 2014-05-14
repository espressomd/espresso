# Non-bonded interactions

cdef class NonBondedInteraction(object):
  
  cdef public object _partTypes
  cdef object _params

  def __init__(self, *args, **kwargs):
    """Represents an instance of a non-bonded interaction, such as lennard jones
    Either called with two particle type id, in which case, the interaction 
    will represent the bonded interaction as it is defined in Espresso core
    Or called with keyword arguments describing a new interaction."""
    
    
    # Interaction id as argument
    if len(args)==2 and isinstance(args[0],int) and isinstance(args[1],int):
      self._partTypes=args
      
      # Load the parameters currently set in the Espresso core
      self._params=self._getParamsFromEsCore()
    
    # Or have we been called with keyword args describing the interaction
    elif len(args)==0:
      # Initialize default values
      self._params=self.defaultParams()
      self._partTypes=[-1,-1]
      
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
      
      self._params = kwargs
      
      # Validation of parameters
      self.validateParams()
      
    else: 
      raise Exception("The constructor has to be called either with two particle type ids (as interger), or with a set of keyword arguments describing a new interaction")

   
      
   


  def isValid(self):
    """Check, if the data stored in the instance still matches what is in Espresso"""

    # check, if the bond parameters saved in the class still match those saved in Espresso
    tempParams =self._getParamsFromEsCore()
    if self._params != tempParams:
      return False
    
    # If we're still here, the instance is valid
    return True
 
 
  def getParams(self):
    """Get interaction parameters"""
    # If this instance refers to an actual interaction defined in the es core, load
    # current parameters from there
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    
    return self._params

  def setParams(self,**p):
    """Update parameters. Only given """
    # Check, if any key was passed, which is not known
    for k in p.keys():
      if k not in self.validKeys():
        raise ValueError("Only the following keys are supported: "+self.validKeys().__str__())

    # When an interaction is newly activated, all required keys must be given
    if not self.isActive():
      for k in self.requiredKeys():
        if k not in p:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
    
    # If this instance refers to an interaction defined in the espresso core,
    # load the parameters from there

    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    
    # Put in values given by the user
    self._params.update(p)
    
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._setParamsInEsCore()

  def validateParams(self):
    return True

  def _getParamsFromEsCore(self):
    raise Exception("Subclasses of NonBondedInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of NonBondedInteraction must define the setParamsFromEsCore() method.")
  
  def defaultParams(self):
    raise Exception("Subclasses of NonBondedInteraction must define the defaultParams() method.")
  
  def isActive(self):
    # If this instance refers to an actual interaction defined in the es core, load
    # current parameters from there
    if self._partTypes[0]>=0 and self._partTypes[1]>=0:
      self._params=self._getParamsFromEsCore()
    raise Exception("Subclasses of NonBondedInteraction must define the isActive() method.")

  
  def typeName(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the requiredKeys() method.")




# Lennard Jones

cdef class LennardJonesInteraction(NonBondedInteraction):
  def validateParams(self):
    if self._params["epsilon"]<0:
      raise ValueError("Lennard-Jones eps has to be >=0")
    if self._params["sigma"]<0:
      raise ValueError("Lennard-Jones sigma has to be >=0")
    if self._params["cutoff"]<0:
      raise ValueError("Lennard-Jones cutoff has to be >=0")
    return True

  def _getParamsFromEsCore(self):
    cdef IA_parameters* iaParams
    iaParams =  get_ia_param(self._partTypes[0],self._partTypes[1]) 
    return {\
      "epsilon":iaParams.LJ_eps,\
      "sigma":iaParams.LJ_sig,\
      "cutoff":iaParams.LJ_cut,\
      "shift":iaParams.LJ_shift,\
      "offset":iaParams.LJ_offset,\
      "min":iaParams.LJ_min}
       

  def isActive(self):
    return (self._params["epsilon"] >0)
  
  def _setParamsInEsCore(self):
    # Handle the case of shift="auto"
    if self._params["shift"]=="auto": 
      # Calc shift
      self._params["shift"]= -( (self._params["sigma"]/self._params["cutoff"])**12 - (self._params["sigma"]/self._params["cutoff"])**6 )
    
    if lennard_jones_set_params(self._partTypes[0],self._partTypes[1],\
                                        self._params["epsilon"], \
                                        self._params["sigma"], \
                                        self._params["cutoff"], \
                                        self._params["shift"], \
                                        self._params["offset"], \
					0.0, \
                                        self._params["min"]):
      raise Exception("Could not set Lennard Jones parameters")					
  
  def defaultParams(self):
    self._params={\
      "epsilon":0.,\
      "sigma":0.,\
      "cutoff":0.,\
      "shift":0.,\
      "offset":0.,\
      "min":0.}

  def typeName(self): 
    return "LennardJones" 
  
  def validKeys(self): 
    return "epsilon","sigma","cutoff","shift","offset","min"
  
  def requiredKeys(self): 
    return "epsilon","sigma","cutoff","shift" 










class NonBondedInteractionHandle(object):
  """Provides access to all Non-bonded interactions between 
  two particle types.""" 

  type1=-1
  type2=-1

  # Here, one line per non-bonded ia
  lennardJones=None


  def __init__(self, _type1, _type2):
    """Takes two particle types as argument"""
    if not (isinstance(_type1,int) and isinstance(_type2,int)):
      raise TypeError("The particle types have to be of type integer.")
    self.type1=_type1
    self.type2=_type2
    
    
    # Here, add one line for each nonbonded ia
    self.lennardJones =LennardJonesInteraction(_type1,_type2)

  
  


cdef class NonBondedInteractions:
  """Access to non-bonded interaction parameters via [i,j], where i,j are particle 
  types. Returns NonBondedInteractionHandle.
  Also: access to force capping
  """
  def __getitem__(self,key):
    if not isinstance(key,tuple):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    if len(key) != 2 or (not isinstance(key[0],int)) or (not isinstance(key[1],int)):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    return NonBondedInteractionHandle(key[0],key[1])
    
  def setForceCap(self,cap):
   if forcecap_set_params(cap):
     raise Exception("Could not set forcecap")

  def getForceCap(self):
    return force_cap






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
      
      self._bondId=bondId
      # Load the parameters currently set in the Espresso core
      self._params=self._getParamsFromEsCore()
      self._bondId=bondId
    
    # Or have we been called with keyword args describing the interaction
    elif len(args)==0:
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
      
      self.params = kwargs
      
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
      for k in p.keys():
        if k not in self.validKeys():
          raise ValueError("Only the following keys are supported: "+self.validKeys().__str__)
      
      # Initialize default values
      self.setDefaultParams()
      # Put in values given by the user
      self._params.update(p)

  def validateParams(self):
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
   fene_set_params(self._bondId,self._params["k"],self._params["d_r_max"],self._params["r_0"])

class HarmonicBond(BondedInteraction):
  def typeNumber(self):
    return 1

  def typeName(self): 
    return "HARMONIC"

  def validKeys(self):
    return "k","r_0","r_cut"

  def requiredKeys(self): 
    return "k","r_0"

  def setDefaultParams(self):
    self._params = {"k'":0.,"r_0":0.,"r_cut":0.} 

  def _getParamsFromEsCore(self):
    return \
      {"k":bonded_ia_params[self._bondId].p.harmonic.k,\
       "r_0":bonded_ia_params[self._bondId].p.harmonic.r,\
       "r_cut":bonded_ia_params[self._bondId].p.harmonic.r_cut}

  def _setParamsInEsCore(self):
   harmonic_set_params(self._bondId,self._params["k"],self._params["r_0"],self._params["r_cut"])

    


bondedInteractionClasses = {0:FeneBond, 1:HarmonicBond}






class BondedInteractions:
  """Represents the non-bonded interactions. Individual interactions can be accessed using
  NonBondedInteractions[i], where i is the bond id. Will return an instance o
  BondedInteractionHandle"""
  def __getitem__(self, key):
    if not isinstance(key,int):
      raise ValueError("Index to BondedInteractions[] hast to ba an integer referring to a bond id")

    # Find out the type of the interaction from Espresso
    bondType = bonded_ia_params[key].type

    # Check if the bonded interaction exists in Espresso core
    if bondType == -1:
      raise ValueError("The bonded interaction with the id "+str(key)+" is not yet defined.")

    # Find the appropriate class representing such a bond
    bondClass =bondedInteractionClasses[bondType]

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






