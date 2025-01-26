#include "integrators/metropolis_hastings_integrator.h"
#include "integrators/monte_carlo_integrator.h"
#include "domain/hypersphere.h"
#include "domain/polygon2d.h"

#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <cmath>

const std::vector<size_t> numPointsValues = {10'000, 100'000, 1'000'000, 10'000'000, 100'000'000, 1'000'000'000};
unsigned int numThreads;

void circleIntegration() {
    //Function to integrate on the circle: x^2 + y^2
    auto f = [](const std::vector<double> &x) {
        return x[0] * x[0] + x[1] * x[1];
    };

    // Circle with radius 1 centered at (0, 0)
    Hypersphere sphere(2, 1.0);

    // Monte Carlo integrator initialization
    MonteCarloIntegrator mcIntegrator(sphere);

    // Print header
    std::cout << "Integrating f(x,y) = x^2 + y^2 over the unit circle (radius = 1):\n";
    std::cout << "Expected result (π/2): " << std::fixed << std::setprecision(6) << (M_PI / 2) << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Time Strat (ms)"
              << std::setw(18) << "Std Result"
              << std::setw(20) << "Strat Result" << "\n";
    std::cout << std::string(95, '-') << "\n";

    //Different sizes for the layered method                                        
    std::vector<int32_t> strataPerDimValues = {5, 10, 20, 50};
    //std::vector<int32_t> strataPerDimValues = {10, 100, 1000, 5000};        //trying something different

    for (size_t numPoints : numPointsValues) {
        // Standard method
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);

        // Printing the standard method line
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"  //No grid for the standard method
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << "-"  // No time for the layered method
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard
                  << std::setw(20) << "-" << "\n";

        // Stratified method
        for (int32_t strataPerDim : strataPerDimValues) {
            auto startStrat = std::chrono::high_resolution_clock::now();
            double resultStratified = mcIntegrator.integrateStratified(f, numPoints, numThreads, strataPerDim);
            auto endStrat = std::chrono::high_resolution_clock::now();
            auto durationStratified = std::chrono::duration_cast<std::chrono::milliseconds>(endStrat - startStrat);

            // Printing of the layered method line
            std::cout << std::setw(12) << numPoints
                      << std::setw(12) << strataPerDim
                      << std::setw(18) << "-"  // No time for the standard method
                      << std::setw(18) << durationStratified.count()
                      << std::setw(18) << "-"
                      << std::setw(20) << resultStratified << "\n";
        }
    }
}

void triangleIntegration() {
    // Integratiion over the (1,1,1) triangle domain with f(x,y)=1
    // to find the area of this triangle.
    // The vertices of an equilateral triangle of side length 1 are as follows:
    // A=(0,0), B=(1,0), C=(0.5, sqrt(3)/2)
    std::vector<std::pair<double, double>> triangleVertices = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.5, std::sqrt(3.0) / 2.0}
    };
    Polygon2D triangle(triangleVertices);

    // Monte Carlo integrator for the triangle
    MonteCarloIntegrator mcIntegrator(triangle);

    // The function to integrate is f(x,y)=1
    auto f = [](const auto &) {
        return 1.0;
    };

    // The expected result is the area of the triangle = sqrt(3)/4
    double expectedTriangleArea = std::sqrt(3.0) / 4.0;

    // Print header
    std::cout << "\nIntegrating f(x,y)=1 over the equilateral triangle (1,1,1):\n";
    std::cout << "Expected area: " << std::fixed << std::setprecision(6) << expectedTriangleArea << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Time Strat (ms)"
              << std::setw(18) << "Std Result"
              << std::setw(20) << "Strat Result" << "\n";
    std::cout << std::string(95, '-') << "\n";

    // Different sizes for the layered method
    std::vector<int32_t> strataPerDimValues = {5, 10, 20, 50};
    //std::vector<int32_t> strataPerDimValues = {10, 100, 1000, 5000};     //trying something different
 
    for (size_t numPoints : numPointsValues) {
        // Standard method
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);

        // Printing the line for the standard method
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"                     // No grid for the standard method
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << "-"                       // No time for the layered method
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard
                  << std::setw(20) << "-" << "\n";

        // Stratified Method
        for (int32_t strataPerDim : strataPerDimValues) {
            auto startStrat = std::chrono::high_resolution_clock::now();
            double resultStratified = mcIntegrator.integrateStratified(f, numPoints, numThreads, strataPerDim);
            auto endStrat = std::chrono::high_resolution_clock::now();
            auto durationStratified = std::chrono::duration_cast<std::chrono::milliseconds>(endStrat - startStrat);

            // Printing the row for the layered method
            std::cout << std::setw(12) << numPoints
                      << std::setw(12) << strataPerDim
                      << std::setw(18) << "-"             // No time for the standard method
                      << std::setw(18) << durationStratified.count()
                      << std::setw(18) << "-"
                      << std::setw(20) << std::fixed << std::setprecision(6) << resultStratified << "\n";
        }
    }
}

void MHintegration() {
    // Suppose we have a standard 2D normal distribution
    auto p = [](const std::vector<double> &x) -> double {
        return 1.0 / (2 * M_PI) * exp(-(x[0] * x[0] + x[1] * x[1]) / 2);
    };

    // Let's calculate E[x_0^2]
    // It should be equal to 1.0 for standard normal distribution
    auto f = [](const std::vector<double> &x) -> double {
        return x[0] * x[0];
    };

    // Integrate the square with the most probability density
    Polygon2D domain({
        {-10.0, -10.0},
        {-10.0, 10.0},
        {10.0, 10.0},
        {10.0, -10.0}
    });

    // MH integrator for the square (set stddev for the Gaussian noise to 0.5)
    MetropolisHastingsIntegrator mh(domain, 0.5);
    std::vector<double> initialPoint = {0.0, 0.0};

    std::cout << "\nCalculating E[x_0^2] of a standard 2D normal distribution:\n";
    std::cout << "Expected result: " << std::fixed << std::setprecision(4) << 1.0 << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
            << std::setw(12) << "Time (ms)"
            << std::setw(13) << "Result"
            << std::setw(18) << "Acceptance rate" << "\n";
    std::cout << std::string(83, '-') << "\n";

    for (size_t numPoints: numPointsValues) {
        // MH is more computationally-heavy
        size_t points = numPoints / 10;

        // Metropolis-Hastings integration on the square
        auto start = std::chrono::high_resolution_clock::now();
        auto [result, acceptanceRate] = mh.integrateParallel(f, p, points, initialPoint, numThreads);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << std::setw(12) << points
                << std::setw(11) << duration.count()
                << std::setw(14) << std::fixed << std::setprecision(6) << result
                << std::setw(18) << std::fixed << std::setprecision(2) << acceptanceRate * 100 << "%\n";
    }
}


int main() {
    numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) {
        // fallback
        numThreads = 16;
    }
    std::cout << "Using " << numThreads << " threads.\n";

    circleIntegration();
    triangleIntegration();
    MHintegration();

    return 0;
}