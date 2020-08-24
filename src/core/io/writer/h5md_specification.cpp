#include "h5md_specification.hpp"
#include "hdf5.h"

namespace Writer {
namespace H5md {

std::array<H5MD_Specification::Dataset, 30> H5MD_Specification::DATASETS = {{
    {"particles/atoms/box/edges", "value", 2, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/box/edges", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/box/edges", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/mass", "value", 2, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/mass", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/mass", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/charge", "value", 2, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/charge", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/charge", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/id", "value", 2, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/id", "step", 1, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/id", "time", 1, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/species", "value", 2, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/species", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/species", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/position", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/position", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/position", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/velocity", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/velocity", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/velocity", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/force", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/force", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/force", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/image", "value", 3, H5T_NATIVE_INT, 3, false},
    {"particles/atoms/image", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/image", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"connectivity/atoms", "value", 3, H5T_NATIVE_INT, 2, false},
    {"connectivity/atoms", "step", 1, H5T_NATIVE_INT, 1, true},
    {"connectivity/atoms", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
}};
}
} // namespace Writer
