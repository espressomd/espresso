#ifndef __VARIANT__HPP
#define __VARIANT__HPP

#include <vector>
#include <string>
#include "Vector.hpp"

class Variant {
public:
  enum Type { NONE = 0, INT, DOUBLE, STRING, INT_VECTOR, DOUBLE_VECTOR };
  Variant() : m_mem(0) { m_type = NONE; }
  Variant(Type t) {
    m_type = t;
    init_mem();
  }
  Variant(const int v) : Variant() {
    construct_from<int>(v);
    m_type = INT;
  }
  Variant(const double v) : Variant() {
    construct_from<double>(v);
    m_type = DOUBLE;
  }
  
  Variant(const std::string &v) : Variant() {
    construct_from<std::string>(v);
    m_type = STRING;
  }

  Variant(std::string &&v) : Variant() {
    construct_from<std::string>(v);
    m_type = STRING;
  }
  
  Variant(const std::vector<int> &v) : Variant() {
    construct_from<typename std::vector<int> >(v);
    m_type = INT_VECTOR;
  }

  Variant(std::vector<int> &&v) : Variant() {
    construct_from<typename std::vector<int> >(v);
    m_type = INT_VECTOR;
  }
  
  Variant(const std::vector<double> &v) : Variant() {
    construct_from<typename std::vector<double> >(v);
    m_type = DOUBLE_VECTOR;
  }

  Variant(std::vector<double> &&v) : Variant() {
    construct_from<typename std::vector<double> >(v);
    m_type = DOUBLE_VECTOR;
  }

  Variant(const Vector3d &v) : Variant() {
    delete_mem();
    std::vector<double> *p = new std::vector<double>;
    m_mem = reinterpret_cast<void *>(p);
    p->resize(3);

    for(int i = 0; i < 3; i++) {
      (*p)[i] = v[i];
    }
    m_type = DOUBLE_VECTOR;
  }
  
  Variant(const Variant &rhs) : Variant() {
    copy(rhs);
  }
  
  ~Variant() {
    delete_mem();
  }
  
  Variant &operator=(const Variant &rhs) {
    copy(rhs);
    return *this;
  }
  
  friend std::ostream &operator<<(std::ostream &os, const Variant &rhs);

  std::uint32_t size() const;
  std::vector<char> serialize() const;
  
  const Type type() const { return m_type; };
  
  operator int()  const {
    if(m_type != INT)
      throw std::string("Variant is not an INT\n");
    return *reinterpret_cast<const int *>(m_mem);
  }
  
  operator double() const {
    if(m_type != DOUBLE)
      throw std::string("Variant is not an DOUBLE\n");
    return *reinterpret_cast<const double *>(m_mem);
  }

  operator std::string&() const {
    if(m_type != STRING)
      throw std::string("Variant is not an STRING\n");
    return *reinterpret_cast<std::string *>(m_mem);
  }

  operator std::vector<int>&() const {
    if(m_type != INT_VECTOR)
      throw std::string("Variant is not an INT_VECTOR\n");
    return *reinterpret_cast<std::vector<int> *>(m_mem);
  }

  operator int*() const {
    if(m_type != INT_VECTOR)
      throw std::string("Variant is not an INT_VECTOR\n");
    return &(*(reinterpret_cast<std::vector<int> *>(m_mem)->begin()));
  }

  operator std::vector<double>&() const {
    if(m_type != DOUBLE_VECTOR)
      throw std::string("Variant is not an DOUBLE_VECTOR\n");
    return *reinterpret_cast<std::vector<double> *>(m_mem);
  }

  operator double*() const {
    if(m_type != DOUBLE_VECTOR)
      throw std::string("Variant is not an DOUBLE_VECTOR\n");
    return &(*(reinterpret_cast<std::vector<double> *>(m_mem)->begin()));
  }

  double &operator[](int i) {
    if(m_type != DOUBLE_VECTOR)
      throw std::string("Variant is not an DOUBLE_VECTOR\n");
    std::vector<double> *d = reinterpret_cast<std::vector<double> *>(m_mem);
    if(i >= d->size())
      d->resize(i+1);
    return (*d)[i];
  }
  
private:
  void copy(const Variant &rhs) {
    switch(rhs.m_type) {
    case NONE:      
      break;
    case INT:
      construct_from<int>(*reinterpret_cast<int *>(rhs.m_mem));      
      break;
    case DOUBLE:
      construct_from<double>(*reinterpret_cast<double *>(rhs.m_mem));
      break;
    case STRING:
      construct_from<std::string>(*reinterpret_cast<std::string *>(rhs.m_mem));
      break;
    case  INT_VECTOR:
      construct_from<typename std::vector<int> >(*reinterpret_cast<typename std::vector<int> *>(rhs.m_mem));
      break;
    case DOUBLE_VECTOR:
      construct_from<typename std::vector<double> >(*reinterpret_cast<typename std::vector<double> *>(rhs.m_mem));
      break;
    }
    m_type = rhs.m_type;
  }
  
  template<typename T>
  void construct_from(const T &v) {
    delete_mem();
    m_mem = reinterpret_cast<void *>(new T(v));
  }

  template<typename T>
  void construct_from(T &&v) {
    delete_mem();
    m_mem = reinterpret_cast<void *>(new T(v));
  }
  
  void init_mem() {
    switch(m_type) {
    case NONE:
      break;
    case INT:
      m_mem = new int;
      break;
    case DOUBLE:
      m_mem = new double;
      break;
    case STRING:
      m_mem = new std::string;
      break;
    case  INT_VECTOR:
      m_mem = new std::vector<int>;
      break;
    case DOUBLE_VECTOR:
      m_mem = new std::vector<double>;
      break;
    }    
  }

  void delete_mem() {
    switch(m_type) {
    case NONE:
      break;
    case INT:
      delete reinterpret_cast<const int *>(m_mem);
      break;
    case DOUBLE:
      delete reinterpret_cast<const double *>(m_mem);
      break;
    case STRING:
      delete reinterpret_cast<const std::string *>(m_mem);
      break;
    case  INT_VECTOR:
      delete reinterpret_cast<const std::vector<int> *>(m_mem);
      break;
    case DOUBLE_VECTOR:
      delete reinterpret_cast<const std::vector<double> *>(m_mem);
      break;
    }
    m_mem = 0;
  }
  Type m_type;
  void *m_mem;
};

#endif
