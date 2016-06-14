#include "TclScriptInterface.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;

/* Source: http://stackoverflow.com/a/6693088/3198615 */
namespace std {
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  std::ostringstream ss;

  std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(ss, " "));
  ss << v.back();

  out << ss.str();
  
  return out;
}
}

namespace ScriptInterface { namespace Tcl {

string TclScriptInterface::print_to_string() const {
  ostringstream res;

  for(auto &p: m_so.get_parameters()) {
    res << p.first << " " << p.second << " ";
  }

  return res.str();
}

void TclScriptInterface::parse_from_string(list<string> &argv) {  
  ParameterMap p = m_so.all_parameters();
  std::map<std::string, Variant> values;
  
  for(list<string>::iterator it = argv.begin(); it != argv.end();) {    
    string s = *it;
    ParameterMap::iterator si = p.find(s);
    if(si == p.end()) {
      ++it;
      continue;
    }
    it = argv.erase(it);
    Parameter p = si->second;

    switch(p.type()) {
      case ParameterType::BOOL:
        values[si->first] = true;
        break;
      case ParameterType::INT:
        {
          stringstream ss(*it);
          int i;
          ss >> i;
          if(ss.fail()) {
            ostringstream error;
            error << s << " expects one integer argument, but got '" << *it << "'";
            throw std::invalid_argument(error.str());
          }
          else {
            values[si->first] = i;
          }
        }
        it = argv.erase(it);
        break;
      case ParameterType::DOUBLE:
        {
          stringstream ss(*it);
          double d;
          ss >> d;
          if(ss.fail()) {
            ostringstream error;
            error << s << " expects one float argument, but got '" << *it << "'";
            throw std::invalid_argument(error.str());
          }
          else {
            values[si->first] = d;
          }
        }
        it = argv.erase(it);
        break;
      case ParameterType::STRING:
        values[si->first] = *it;
        it = argv.erase(it);
        break;
      case ParameterType::INT_VECTOR:
        {
          vector<int> v;
          for(int j = 1; j <= p.n_elements(); ++j) {
            {
              stringstream ss(*it);
              int i;
              ss >> i;
              if(ss.fail()) {
                ostringstream error;
                error << s << " expects " << p.n_elements() << " integer arguments, but argument " << j << " was '" << *it << "'";
                throw std::invalid_argument(error.str());
              }
              else
                v.push_back(i);

              it = argv.erase(it);
            }
          }
          values[si->first] = v;
          break;
        }
      case ParameterType::DOUBLE_VECTOR:
        {
          vector<double> v;
          for(int j = 1; j <= p.n_elements(); ++j) {
            {
              stringstream ss(*it);
              double i;
              ss >> i;
              if(ss.fail()) {
                ostringstream error;
                error << s << " expects " << p.n_elements() << " float arguments, but argument " << j << " was '" << *it << "'";
                throw std::invalid_argument(error.str());
              }
              else
                v.push_back(i);

              it = argv.erase(it);
            }
          }
          values[si->first] = v;
          break;
        }
    }
  }

  /* Forward the parameters to m_so */
  m_so.set_parameters(values);
}

}} /* namespace */

