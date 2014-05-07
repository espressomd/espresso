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
      self._partTypes=[-1,-1]
      # Check if all required keys are given
      for k in self.requiredKeys():
        if k not in kwargs:
          raise ValueError("At least the following keys have to be given as keyword arguments: "+self.requiredKeys().__str__())
      
      self.params = kwargs
      
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
    raise Exception("Subclasses of NonBondedInteraction must define the _getParamsFromEsCore() method.")
  
  def _setParamsInEsCore(self):
    raise Exception("Subclasses of NonBondedInteraction must define the setParamsFromEsCore() method.")
  
  def setDefaultParams(self):
    raise Exception("Subclasses of NonBondedInteraction must define the setDefaultParams() method.")

  
  def typeName(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the typeName() method.")
  
  def validKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the validKeys() method.")
  
  def requiredKeys(self): 
    raise Exception("Subclasses of NonBondedInteraction must define the requiredKeys() method.")




# Lennard Jones

cdef class LennardJonesInteraction(NonBondedInteraction):
  def validateParams(self):
    if self.params["epsilon"]<0:
      raise ValueError("Lennard-Jones eps has to be >=0")
    if self.params["sigma"]<0:
      raise ValueError("Lennard-Jones sigma has to be >=0")
    if self.params["cutoff"]<0:
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
      "capradius":iaParams.LJ_capradius}
       

    
  
  def _setParamsInEsCore(self):
    if lennard_jones_set_params(self._partTypes[0],self._partTypes[1],\
                                        self._params["epsilon"], \
                                        self._params["sigma"], \
                                        self._params["cutoff"], \
                                        self._params["shift"], \
                                        self._params["offset"], \
                                        self._params["capradius"], \
                                        self._params["min"]):
      raise Exception("Could not set Lennard Jones parameters")					
    ljforcecap_set_params(self._params["capradius"])
  
  def setDefaultParams(self):
    self._params={\
      "epsilon":0.,\
      "sigma":0.,\
      "cutoff":0.,\
      "shift":0.,\
      "offset":0.,\
      "capradius":0.,\
      "min":0.}

  def typeName(self): 
    return "LennardJones" 
  
  def validKeys(self): 
    return "epsilon","sigma","cutoff","shift","offset","capradius","min"
  
  def requiredKeys(self): 
    return "epsilon","sigma","cutoff","shift" 










class NonBondedInteractionHandle(object):
  """Provides access to all Non-bonded interactions between 
  two particle types.""" 

  type1=-1
  type2=-1

  def __init__(self, _type1, _type2):
    """Takes two particle types as argument"""
    if not (isinstance(_type1,int) and isinstance(_type2,int)):
      raise TypeError("The particle types have to be of type integer.")
    self.type1=_type1
    self.type2=_type2
    
    # Here, add one line for each nonbonded ia
    self.__class__.lennardJones =self.iaProperty(LennardJonesInteraction,"LennardJones")

  def _getter(self,iaClass):
    """Generates a getter funciton for a speciffic interaction type such as Lennard Jones""" 
    return lambda self: iaClass(self.type1,self.type2)

  def _setter(self,_iaClass):
    """Generates a setter funciton for a speciffic interaction type such as Lennard Jones""" 

    iaClass = _iaClass

    def f(self,v):
      if not isinstance(v,iaClass):
        raise TypeError("The value given must be of type "+iaClass.__name__)
      v._partTypes[0]=self.type1
      v._partTypes[1]=self.type2
      v._setParamsInEsCore()

    return f

  def iaProperty(self,iaClass,iaName):
    """Generate a property with the correct setter and getter for
    a specific interaction such as Lennard Jones"""
    return property(fget=self._getter(iaClass),fset=self._setter(iaClass),doc=iaName)
  
  


cdef class NonBondedInteractions:
  """Access to non-bonded interaction parameters via [i,j], where i,j are particle 
  types. Returns NonBondedInteractionHandle."""
  def __getitem__(self,key):
    if not isinstance(key,tuple):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    if len(key) != 2 or (not isinstance(key[0],int)) or (not isinstance(key[1],int)):
      raise ValueError("NonBondedInteractions[] expects two particle types as indices.")
    return NonBondedInteractionHandle(key[0],key[1])







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






