#ifndef METROPOLIS_HASTINGS_ISING_H
#define METROPOLIS_HASTINGS_ISING_H

#include "abstract_integrator.h"
#include "../domain/dummy_domain.h"

#include <random>
#include <functional>
#include <vector>
#include <future>
#include <cstdint>

// Metropolis-Hastings for ising model
class MetropolisHastingsIsing : public AbstractIntegrator {
public:
    // Constructor
    explicit MetropolisHastingsIsing();

    // Perform a single Metropolis-Hastings chain for the Ising model
    // Returns a pair of {energy average, #accepted points}
    std::pair<double, int32_t> integrateSingleChainIsing(
        const std::function<double(const std::vector<std::vector<int>> &)> &f,
        size_t numPoints,
        std::vector<std::vector<int>> &lattice,
        double temperature,
        std::mt19937 &engine
    );

    // Perform parallel computation of energy average for the Ising model using multiple chains
    // Returns a pair of {energy average, acceptance rate}
    std::pair<double, double> integrateParallelIsing(
        const std::function<double(const std::vector<std::vector<int>> &)> &f,
        size_t numPoints,
        std::vector<std::vector<int>> &initialLattice,
        double temperature,
        size_t numChains
    );

private:
    double calculateEnergyDifference(const std::vector<std::vector<int>> &lattice, int32_t row, int32_t col);
    void flipSpin(std::vector<std::vector<int>> &lattice, int32_t &row, int32_t &col, std::mt19937 &engine);
};

#endif // METROPOLIS_HASTINGS_ISING_H