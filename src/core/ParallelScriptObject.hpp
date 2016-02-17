#ifndef __PARALLEL_SCRIPTOBJECTP_HPP
#define __PARALLEL_SCRIPTOBJECTP_HPP

#include "ScriptObject.hpp"
#include "ParallelObject.hpp"

class ParallelScriptObject : public ScriptObject, public ParallelObject {
 public:
  virtual void set_parameters_all_nodes(Parameters &parameters);
  virtual ~ParallelScriptObject();
  
 private:
  virtual void callback(int par);   
  enum CallbackActions { DELETE, SET_PARAMETERS, SET_PARAMETER };  
};

#endif
