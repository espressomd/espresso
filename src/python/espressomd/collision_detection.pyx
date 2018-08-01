from .script_interface import ScriptInterfaceHelper,script_interface_register
from .utils import handle_errors,to_str
from .interactions import BondedInteraction,BondedInteractions


cdef extern from "collision.hpp":
    const int COLLISION_MODE_OFF
    const int COLLISION_MODE_BOND
    const int COLLISION_MODE_VS 
    const int COLLISION_MODE_GLUE_TO_SURF  
    const int COLLISION_MODE_BIND_THREE_PARTICLES  
    

@script_interface_register
class CollisionDetection(ScriptInterfaceHelper):
    """Inteface to the collision detection / dynamic binding.
    
       See :ref:`Creating bonds when particles collide` for detailed instructions.
    
    """

    _so_name = "CollisionDetection::CollisionDetection"

    def __init__(self,*args,**kwargs):
        """
        This class should not be instanced by the user. Instead, use
        the :attr:`espressomd.system.System.collision_detection` attribute
        of the system class to access the collision detection.
           
        Use :meth:`espressomd.collision_detection.CollisionDetection.set_params`
        to change the parameters of the collision detection.

        """
        # If no mode is specified at construction, use off.
        if "mode" not in kwargs:
            kwargs["mode"]="off"
        super(type(self),self).__init__()
        self.set_params(**kwargs)
        

    def validate(self):
        """Validates the parameters of the collision detection.

           This is called automatically on parameter change

        """
        return self.call_method("validate")
    
    # Do not allow setting of individual attributes
    def __setattr__(self,*args,**kwargs):
        raise Exception("Please et all parameters at once via collision_detection.set_params()")

    # Override to call validat after parameter update
    def set_params(self, **kwargs):
        """
        Set the parameters for the collision detection
           
        See :ref:`Creating bonds when particles collide` for detailed instructions.


        Parameters
        ----------
        mode : One of "off", "bind_centers", "bind_at_point_of_collision", "bind_three_particles", "glue_to_surface"
               Collision deteciton mode
          
        distance : :obj:`float`
               Distance below which a pair of particles is considered in the collision detection
          
        bond_centers : Instance of :class:`espressomd.interactions.BondedInteraction`
               Bond to add between the colliding particles
          
        bond_vs :  Instance of :class:`espressomd.interactions.BondedInteraction`
               Bond to add between virtual sites (for modes using virtual sites)
           
        part_type_vs : :obj:`int`
               Particle type of the virtual sites being created on collision (virtual sites based modes)
           
        part_type_to_be_glued : :obj:`int`
               particle type for "glue_to_surface|" mode. See user guide.
           
        part_type_to_attach_vs_to : :obj:`int`
               particle type for "glue_to_surface|" mode. See user guide.
           
        part_type_after_glueing : :obj:`int`
               particle type for "glue_to_surface|" mode. See user guide.
           
        distance_glued_particle_to_vs : :obj:`float`
               Distnace for "glue_to_surface" mode. See user guide.
           
        bond_three_particles : Instance of :class:`espressomd.interactions.BondedInteraction`
               First angular bond for the "bind_three_particles" mode. See user guide
          
        three_particle_binding_angle_resolution : :obj:`int`
              Resolution for the angular bonds (mode "bind_three_particles"). 
              Resolution+1 bonds are needed to accomodate the case of a 180 degrees

        """


        if not ("mode" in kwargs):
            raise Exception("Collision mode must be specified via the mode keyword argument")
        
        # Completeness of parameter set
        if not (set(kwargs.keys()) == set(self._params_for_mode(kwargs["mode"]))):
            raise Exception("Parameter set does not match mode. ",kwargs["mode"],"requries ",self._params_for_mode(kwargs["mode"]))

        # Mode
        kwargs["mode"]=self._int_mode[kwargs["mode"]]

        # Convert bonds to bond ids
        for name in [ "bond_centers","bond_vs","bond_three_particle_binding"]:
            if name in kwargs:
                if isinstance(kwargs[name],BondedInteraction):
                    kwargs[name]=kwargs[name]._bond_id
        super(type(self),self).set_params(**kwargs)
        self.validate()
        handle_errors("Validation of collision detection failed")

    def get_parameter(self,name):
        #"""Gets a single parameter from the collision detection."""
        
        res=super(type(self),self).get_parameter(name)
        return self._convert_param(name,res)
    
    def get_params(self):
        """Returns the parameters of the collision detection as dict.
        
        """
        res=super(type(self),self).get_params()
        for k in res.keys():
            res[k]=self._convert_param(k,res[k])
        
        # Filter key-value paris according to active mode
        return {k:res[k] for k in self._params_for_mode(res["mode"])}


    def _convert_param(self,name,value):
        """Handles type conversion core -> python
            
            Bond types: int -> BondedInteraction
            mode: int -> string

            """
        # Py3: Cast from binary to normal string. Don't understand, why a
        # binary string can even occur, here, but it does.
        name=to_str(name)
        # Convert int mode parameter to string
        res=value
        if name == "mode":
            res=self._str_mode(value)
        
        # Convert bond parameters from bond ids to into BondedInteractions
        if name in [ "bond_centers","bond_vs","bond_three_particle_binding"]:
            if value==-1: # Not defined
                res=None
            else:
                res=BondedInteractions()[value]
        return res
    
    def _params_for_mode(self,mode):
        """The parameter names expected for a given collision mode
        
        """
        if mode == "off":
            return ("mode",)
        if mode == "bind_centers":
            return ("mode","bond_centers","distance")
        if mode == "bind_at_point_of_collision": 
            return ("mode","bond_centers","bond_vs","part_type_vs","distance","vs_placement")
        if mode == "glue_to_surface":
            return ("mode","bond_centers","bond_vs","part_type_vs","part_type_to_be_glued","part_type_to_attach_vs_to","part_type_after_glueing","distance","distance_glued_particle_to_vs")
        if mode == "bind_three_particles":
            return ("mode","bond_centers","distance","bond_three_particles","three_particle_binding_angle_resolution")
        raise Exception("Mode not hanled: "+mode.__str__())
            
            
    
    _int_mode={
        "off":int(COLLISION_MODE_OFF),
        "bind_centers":int(COLLISION_MODE_BOND),
        "bind_at_point_of_collision":int(COLLISION_MODE_VS),
        "glue_to_surface":int(COLLISION_MODE_GLUE_TO_SURF),
        "bind_three_particles":int(COLLISION_MODE_BIND_THREE_PARTICLES)}
    
    def _str_mode(self,int_mode):
        """String mode name from int ones provided by the core
        
        """
        for key in self._int_mode:
            if self._int_mode[key] == int_mode:
                return key
        raise Exception("Unknown integer collision mode %d" % int_mode)
    
    
