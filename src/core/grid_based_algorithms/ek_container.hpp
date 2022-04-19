#ifndef ESPRESSO_EK_CONTAINER_HPP
#define ESPRESSO_EK_CONTAINER_HPP

#include "EKContainer.hpp"
#include "electrokinetics/EKinWalberlaBase.hpp"

namespace EK {
extern EKContainer<EKinWalberlaBase<double>> ek_container;

double get_tau();
int get_steps_per_md_step(double md_timestep);
void propagate();
} // namespace EK

#endif // ESPRESSO_EK_CONTAINER_HPP
