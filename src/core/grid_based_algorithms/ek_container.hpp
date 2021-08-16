#ifndef ESPRESSO_EK_CONTAINER_HPP
#define ESPRESSO_EK_CONTAINER_HPP

#include "EKContainer.hpp"
#include "EKinWalberlaBase.hpp"

namespace EK {
extern EKContainer<EKinWalberlaBase<double>> ek_container;

double get_tau();
void propagate();
} // namespace EK

#endif // ESPRESSO_EK_CONTAINER_HPP
