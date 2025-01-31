#include "ising_model.h"

// compute total energy of the system
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

// execute ising model and compute specific heat per particle
void isingModel() {
    size_t latticeSize = 20;     // size of the square lattice
    double T = 4.5;              // temperature of the system
    double dT = 0.1;            // increment in T
    size_t numPoints = 1000000;   // points for MCMC
    
    size_t numChains = std::thread::hardware_concurrency();
    if (numChains == 0) {
        // fallback
        numChains = 16;
    }
    std::cout << "Using " << numChains << " threads and " << numPoints << " points.\n";

    // lattice initialization
    std::mt19937 engine;
    engine.seed(std::random_device{}());

    std::vector<std::vector<int>> lattice(latticeSize, std::vector<int>(latticeSize));
    std::uniform_int_distribution<> spinDist(-1, 1);
    for (auto &row : lattice) {
        for (auto &spin : row) {
            spin = (spinDist(engine) == 0) ? 1 : -1;
        }
    }

    MetropolisHastingsIsing mhIsing;

    // anonymous function for improved flexibility and modularity in case we want to change the energy function
    auto energyFunction = [](const std::vector<std::vector<int>> &lat) {
        return calculateTotalEnergy(lat);
    };

    // avg Energy and rate of change with respect to T
    auto [averageEnergy_T, acceptanceRate_T] = mhIsing.integrateParallelIsing(energyFunction, numPoints, lattice, T, numChains);
    auto [averageEnergy_TdT, acceptanceRate_TdT] = mhIsing.integrateParallelIsing(energyFunction, numPoints, lattice, T + dT, numChains);

    // compute rate of change for specific heat per particle
    double specificHeat = (averageEnergy_TdT - averageEnergy_T) / dT;

    // output
    std::cout << "Temperature: " << T << std::endl;
    std::cout << "Average Energy at T: " << averageEnergy_T << std::endl;
    std::cout << "Average Energy at T + dT: " << averageEnergy_TdT << std::endl;
    std::cout << "Specific Heat: " << specificHeat << std::endl;
}