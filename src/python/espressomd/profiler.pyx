cimport c_profiler

def begin_section(name):
    c_profiler.begin_section(name)

def end_section(name):
    c_profiler.end_section(name)
