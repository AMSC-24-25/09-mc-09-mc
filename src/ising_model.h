#ifndef ISING_MODEL_H
#define ISING_MODEL_H

#include "integrators/metropolis_hastings_ising.h"

#include <iostream>
#include <random>
#include <vector>
#include <cstdint>

// total energy of the system
double calculateTotalEnergy(const std::vector<std::vector<int>> &lattice);

// execute ising model and compute specific heat per particle
void isingModel();

#endif