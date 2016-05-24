from actors cimport *
from libcpp.string cimport string  # import std::string

include "myconfig.pxi"


# As cython does not have multiple inheritance for cdef classes
# we don't put the stuff in a class and import it method by method


IF SCAFACOS==1:
    
    def valid_keys(self):
        return "method_name","method_params"
    
    def required_keys(self):
        return "method_name","method_params"
      
    def validate_params(self):
        return True
    
    def _get_params_from_es_core(self):
        # Parameters are returned as strings
        # First word is method name, rest are key value pairs
        # which we convert to a dict
        p=get_parameters().split(" ")
        res={}
        res["method_name"]=p[0]
        for i in range((len(p)-1)/2):
          res[p[2*i+1]]=p[2*i+2]
        return res
    
    def _set_params_in_es_core(self):
        # Convert the parameter dictionary to a list of strings
        method_parameters=""
        for k in self._params["method_params"].keys():
          method_parameters+=k+" "+str(self._params[k])
    
        set_parameters(self._params["method_name"], method_parameters,self.dipolar)
    
    def _activate_method(self):
        raise Exception("Needs to be dfined in subclass.")
    
    def _deactivate_method(self):
        raise Exception("Needs to be dfined in subclass.")
        
        
        
    
    
    
