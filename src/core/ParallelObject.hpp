#ifndef __PARALLELOBJECT_HPP
#define __PARALLELOBJECT_HPP

#include <unordered_map>
#include <memory>

class ParallelObject {
 public:
  ParallelObject();
  int get_id() const;

  /** This probably can and should be a pointer to const */
  static ParallelObject *get_local_address(int id);
 protected:
  void call_slaves(int par);  
 private:
  virtual void callback(int par) = 0;
  int m_callback_id;
  
  static std::unordered_map<int, ParallelObject *> address_table;
};

#endif
