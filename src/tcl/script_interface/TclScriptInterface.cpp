#include "TclScriptInterface.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;

namespace Utils {
std::ostream &
operator<<(std::ostream &out,
           ObjectId<ScriptInterface::ScriptInterfaceBase> const &oid) {
  out << oid.m_id;

  return out;
}
}

namespace ScriptInterface {
namespace Tcl {

template <typename T>
T read_element(std::list<std::string>::iterator &it,
               std::list<std::string> &argv) {
  stringstream ss(*it);
  T value;
  ss >> value;

  if (ss.fail()) {
    throw std::exception();
  }

  it = argv.erase(it);

  return value;
}

template <typename T>
std::vector<T> read_n_elements(int n_elements,
                               std::list<std::string>::iterator &it,
                               std::list<std::string> &argv) {
  std::vector<T> v;

  for (int i = 0; i < n_elements; i++) {
    v.push_back(read_element<T>(it, argv));
  }

  return v;
}

string TclScriptInterface::print_to_string() const {
  ostringstream res;

  for (auto &p : m_so->get_parameters()) {
    res << p.first << " " << p.second << " ";
  }

  return res.str();
}

void TclScriptInterface::parse_from_string(list<string> &argv) {
  ParameterMap p = m_so->valid_parameters();
  std::map<std::string, Variant> values;

  for (list<string>::iterator it = argv.begin(); it != argv.end();) {
    string s = *it;
    ParameterMap::iterator si = p.find(s);
    if (si == p.end()) {
      ++it;
      continue;
    }
    it = argv.erase(it);
    Parameter p = si->second;

    switch (p.type()) {
    case ParameterType::BOOL:
      values[si->first] = true;
      break;
    case ParameterType::INT: {
      try {
        values[si->first] = read_element<int>(it, argv);
      } catch (std::exception &e) {
        throw std::invalid_argument(
            s + " expects one integer argument, but got '" + *it + "'");
      }
      break;
    case ParameterType::DOUBLE:
      try {
        values[si->first] = read_element<double>(it, argv);
      } catch (std::exception &e) {
        throw std::invalid_argument(
            s + " expects one float argument, but got '" + *it + "'");
      }
      break;
    case ParameterType::STRING:
      values[si->first] = *it;
      it = argv.erase(it);
      break;
    case ParameterType::INT_VECTOR:
      try {
        values[si->first] = read_n_elements<int>(p.n_elements(), it, argv);
      } catch (std::exception &e) {
        ostringstream error;
        error << s << " expects " << p.n_elements()
              << " integer arguments, but got '" << *it << "'";
        throw std::invalid_argument(error.str());
      }
      break;

    case ParameterType::VECTOR2D:
      try {
        values[si->first] =
            Vector2d(read_n_elements<double>(p.n_elements(), it, argv));
      } catch (std::exception &e) {
        ostringstream error;
        error << s << " expects " << p.n_elements()
              << " float arguments, but got '" << *it << "'";
        throw std::invalid_argument(error.str());
      }
      break;

    case ParameterType::VECTOR3D:
      try {
        values[si->first] =
            Vector3d(read_n_elements<double>(p.n_elements(), it, argv));
      } catch (std::exception &e) {
        ostringstream error;
        error << s << " expects " << p.n_elements()
              << " float arguments, but got '" << *it << "'";
        throw std::invalid_argument(error.str());
      }
      break;

    case ParameterType::DOUBLE_VECTOR:
      try {
        values[si->first] = read_n_elements<double>(p.n_elements(), it, argv);
      } catch (std::exception &e) {
        ostringstream error;
        error << s << " expects " << p.n_elements()
              << " float arguments, but got '" << *it << "'";
        throw std::invalid_argument(error.str());
      }
      break;
    }
    }
  }

  /* Forward the parameters to m_so */
  m_so->set_parameters(values);
}
}
} /* namespace */
