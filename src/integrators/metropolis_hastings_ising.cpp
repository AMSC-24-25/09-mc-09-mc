#include "metropolis_hastings_ising.h"

// initialize integrator with a dummy domain (lattice is used as domain instead)
MetropolisHastingsIsing::MetropolisHastingsIsing()
    : AbstractIntegrator(DummyDomain()) {}

// flip randomly selected spin to generate new state
void MetropolisHastingsIsing::flipSpin(std::vector<std::vector<int>> &lattice, int32_t &row, int32_t &col, std::mt19937 &engine) {
    std::uniform_int_distribution<> rowDist(0, lattice.size() - 1);
    std::uniform_int_distribution<> colDist(0, lattice[0].size() - 1);
    
    row = rowDist(engine);
    col = colDist(engine);

    lattice[row][col] = -lattice[row][col];
}

// energy difference for acceptance ratio (to avoid computing the energy function twice)
double MetropolisHastingsIsing::calculateEnergyDifference(const std::vector<std::vector<int>> &lattice, int32_t row, int32_t col) {
    size_t rows = lattice.size();
    size_t cols = lattice[0].size();
    int32_t spin = lattice[row][col];
    int32_t rightNeighbor = lattice[row][(col + 1) % cols];
    int32_t leftNeighbor = lattice[row][(col + cols - 1) % cols];
    int32_t topNeighbor = lattice[(row + rows - 1) % rows][col];
    int32_t bottomNeighbor = lattice[(row + 1) % rows][col];

    double deltaE = 2 * spin * (rightNeighbor + leftNeighbor + topNeighbor + bottomNeighbor);
    return deltaE;
}

// perform a single Metropolis-Hastings chain for the Ising model
std::pair<double, int32_t> MetropolisHastingsIsing::integrateSingleChainIsing(
    const std::function<double(const std::vector<std::vector<int>> &)> &f,
    size_t numPoints,
    std::vector<std::vector<int>> &lattice,
    double temperature,
    std::mt19937 &engine
) {
    double sumF = 0.0;
    int32_t accepted = 0;
    int32_t samples = 0;

    std::uniform_real_distribution<> unifDist(0.0, 1.0);

    // main sampling loop
    for (size_t i = 0; i < numPoints; ++i) {
        std::vector<std::vector<int>> candidateLattice = lattice;
        int32_t row, col;
        flipSpin(candidateLattice, row, col, engine);

        // acceptance ratio uses energy difference for probability of changing state
        double deltaE = calculateEnergyDifference(lattice, row, col);
        double acceptanceRatio = exp(-deltaE / temperature);

        if (unifDist(engine) <= acceptanceRatio) {
            lattice = candidateLattice;
            ++accepted;
        }

        sumF += f(lattice);
        ++samples;
    }

    return {sumF / samples, accepted};
}

// perform parallel computation of energy average for the Ising model using multiple chains
std::pair<double, double> MetropolisHastingsIsing::integrateParallelIsing(
    const std::function<double(const std::vector<std::vector<int>> &)> &f,
    size_t numPoints,
    std::vector<std::vector<int>> &initialLattice,
    double temperature,
    size_t numChains
) {
    // initialize random engines through abstractIntegrator
    initializeEngines(numChains);
    size_t pointsPerChain = numPoints / numChains;

    std::vector<std::future<std::pair<double, int32_t>>> futures;
    futures.reserve(numChains);
    for (size_t i = 0; i < numChains; ++i) {
        futures.push_back(std::async(
            std::launch::async,
            [this, &f, pointsPerChain, &initialLattice, temperature, &engine = engines[i]]() {
                std::vector<std::vector<int>> latticeCopy = initialLattice;
                return this->integrateSingleChainIsing(f, pointsPerChain, latticeCopy, temperature, engine);
            }
        ));
    }

    double totalResult = 0.0;
    int32_t totalAccepted = 0;
    for (auto &future: futures) {
        auto [result, accepted] = future.get();
        totalResult += result;
        totalAccepted += accepted;
    }

    double acceptanceRate = static_cast<double>(totalAccepted) / numPoints;

    return {totalResult / numChains, acceptanceRate};
}
