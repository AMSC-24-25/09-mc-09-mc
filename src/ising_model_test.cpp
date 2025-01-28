#include <iostream>
#include <functional>
#include <random>
#include <vector>
#include <cmath>
#include "integrators/metropolis_hastings_ising.h"

// Total energy of the system
double calculateTotalEnergy(const std::vector<std::vector<int>> &lattice) {
    double energy = 0.0;
    for (size_t i = 0; i < lattice.size(); ++i) {
        for (size_t j = 0; j < lattice[i].size(); ++j) {
            int32_t spin = lattice[i][j];
            int32_t rightNeighbor = lattice[i][(j + 1) % lattice[i].size()];
            int32_t bottomNeighbor = lattice[(i + 1) % lattice.size()][j];
            energy -= spin * (rightNeighbor + bottomNeighbor);
        }
    }
    return energy;
}

// Specific heat
void performIsingSimulation() {
    // Simulation parameters
    size_t latticeSize = 10;     // Lattice
    double T = 4.5;              // Temperature
    double dT = 0.1;            // Increment to calculate derivative
    size_t numPoints = 1000000;   // Points for MCMC
    size_t numChains = 16;        // Number of parallel chains

    std::mt19937 engine;
    engine.seed(std::random_device{}());

    // Lattice initialization
    std::vector<std::vector<int>> lattice(latticeSize, std::vector<int>(latticeSize));
    std::uniform_int_distribution<> spinDist(-1, 1);
    for (auto &row : lattice) {
        for (auto &spin : row) {
            spin = (spinDist(engine) == 0) ? 1 : -1;
        }
    }

    MetropolisHastingsIsing mhIsing;

    auto energyFunction = [](const std::vector<std::vector<int>> &lat) {
        return calculateTotalEnergy(lat);
    };

    // Avg Energy change with respect to temperature (specific heat in T)
    auto [averageEnergy_T, acceptanceRate_T] = mhIsing.integrateParallelIsing(energyFunction, numPoints, lattice, T, numChains);
    auto [averageEnergy_TdT, acceptanceRate_TdT] = mhIsing.integrateParallelIsing(energyFunction, numPoints, lattice, T + dT, numChains);

    double specificHeat = (averageEnergy_TdT - averageEnergy_T) / dT;

    // Output
    std::cout << "Temperature: " << T << std::endl;
    std::cout << "Average Energy at T: " << averageEnergy_T << std::endl;
    std::cout << "Average Energy at T + dT: " << averageEnergy_TdT << std::endl;
    std::cout << "Specific Heat: " << specificHeat << std::endl;
}