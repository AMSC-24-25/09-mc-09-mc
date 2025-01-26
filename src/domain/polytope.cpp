#include "polytope.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <numeric>
#include "polytope.h"
#include <stdexcept>

// Constructor
Polytopes::Polytopes(const std::vector<std::vector<double>>& A, const std::vector<double>& b)
    : A(A), b(b) {
    if (A.empty() || b.empty()) {
        throw std::invalid_argument("Matrix A and vector b cannot be empty.");
    }
    if (A.size() != b.size()) {
        throw std::invalid_argument("The number of rows in A must match the size of b.");
    }
} 

// Check if a point lies within the polytope
bool Polytopes::contains(const std::vector<double> &point) const {
    if (point.size() != A[0].size()) {
        throw std::invalid_argument("Point dimensionality does not match the polytope.");
    }

   for (size_t i = 0; i < A.size(); ++i) {
    double dotProduct = std::inner_product(A[i].begin(), A[i].end(), point.begin(), 0.0);

    std::cout << "Testing inequality " << i << ": "
              << "dotProduct = " << dotProduct << ", b[" << i << "] = " << b[i] << "\n";

    if (dotProduct > b[i] + 1e-9) {
        std::cout << "Point (" << point[0] << ", " << point[1] << ") fails inequality " << i << "\n";
        return false;
}
    }
    return true;
}


// Return the dimensions of the polytope
size_t Polytopes::getDimensions() const {
    return A[0].size();
}


std::vector<std::pair<double, double>> Polytopes::getBounds() const {
    size_t dimensions = A[0].size();
    std::vector<std::pair<double, double>> bounds(dimensions, {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()});

    for (size_t dim = 0; dim < dimensions; ++dim) {
        double minBound = std::numeric_limits<double>::lowest();
        double maxBound = std::numeric_limits<double>::max();

        for (size_t i = 0; i < A.size(); ++i) {
            if (A[i][dim] != 0) {  // Only consider constraints involving this dimension
                double bound = b[i] / A[i][dim];

                if (A[i][dim] > 0) {
                    maxBound = std::min(maxBound, bound);
                } else {
                    minBound = std::max(minBound, bound);
                }
            }
        }

        // Validate bounds for this dimension
        if (minBound > maxBound) {
            throw std::runtime_error("Invalid bounds: check the polytope definition.");
        }

        bounds[dim] = {std::max(0.0, minBound), std::min(1.0, maxBound)};
    }

    return bounds;
}






// Return the approximate volume of the polytope (default: -1 for no calculation)
double Polytopes::getBoundedVolume() const {
    return -1.0; // Not used with Metropolis-Hastings
}
