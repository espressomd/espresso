#include <profiler/profiler.hpp>

#include <string>

namespace Profiler {
    void begin_section(const std::string &name) {
        ESPRESSO_PROFILER_MARK_BEGIN(name.c_str());
    }
    void end_section(const std::string &name) {
        ESPRESSO_PROFILER_MARK_END(name.c_str());
    }
}


