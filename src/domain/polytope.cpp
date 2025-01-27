#include "polytope.h"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <cdd/setoper.h>
#include <cdd/cdd.h>
#include <numeric>

// Constructor
Polytopes::Polytopes(const std::vector<std::vector<double>>& vertices)
    : vertices(vertices) {
    if (vertices.empty()) {
        throw std::invalid_argument("The vertices cannot be empty.");
    }
}



// Check if a point lies within the polytope
bool Polytopes::contains(const std::vector<double>& point) const {
    size_t n_points = vertices.size(); // Number of vertices
    size_t n_dim = point.size();       // Dimensionality of the point

    // Ensure the point and vertices have compatible dimensions
    if (n_dim != vertices[0].size()) {
        throw std::invalid_argument("Point dimensionality does not match the polytope.");
    }

    // Formulate the A_eq matrix (transpose of vertices + 1 row of ones)
    Eigen::MatrixXd A_eq(n_dim + 1, n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < n_dim; ++j) {
            A_eq(j, i) = vertices[i][j]; // Transpose of vertices
        }
        A_eq(n_dim, i) = 1.0; // Last row is all ones
    }

    // Formulate the b_eq vector (the target point + 1 for convex combination)
    Eigen::VectorXd b_eq(n_dim + 1);
    for (size_t i = 0; i < n_dim; ++i) {
        b_eq(i) = point[i];
    }
    b_eq(n_dim) = 1.0;

    // Objective function (minimize a zero vector since we're only checking feasibility)
    Eigen::VectorXd c = Eigen::VectorXd::Zero(n_points);

    // Solve the linear programming problem using a simplex algorithm
    Eigen::VectorXd solution(n_points);
    try {
        // Use Eigen's full-pivot solver to find a solution
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(A_eq.transpose());
        solution = solver.solve(b_eq);

        // Check if all coefficients are non-negative (convex combination constraint)
        for (size_t i = 0; i < n_points; ++i) {
            if (solution(i) < -1e-9) {
                return false; // Not a valid convex combination
            }
        }
        return true; // Point is inside the convex hull
    } catch (...) {
        return false; // If any error occurs during the solving process
    }
}

size_t Polytopes::getDimensions() const {
    return vertices.empty() ? 0 : vertices[0].size();
}


std::vector<std::pair<double, double>> Polytopes::getBounds() const {
    size_t dimensions = getDimensions();

    // Initialize bounds for each dimension
    std::vector<std::pair<double, double>> bounds(dimensions, {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()});

    // Iterate through the vertices to compute the min and max for each dimension
    for (const auto& vertex : vertices) {
        for (size_t dim = 0; dim < dimensions; ++dim) {
            bounds[dim].first = std::min(bounds[dim].first, vertex[dim]);  // Update minimum
            bounds[dim].second = std::max(bounds[dim].second, vertex[dim]); // Update maximum
        }
    }

    return bounds;
}





double Polytopes::getBoundedVolume() const {
    // Decompose the polytope into simplexes
    auto simplexes = decomposePolytope(vertices);

    // Calculate the total volume by summing up simplex volumes
    double totalVolume = 0.0;
    for (const auto& simplex : simplexes) {
        double simplexVolume = calculateSimplexVolume(simplex);
        totalVolume += simplexVolume;
    }

    return totalVolume;
}


// Function to decompose a 15D polytope into simplexes
std::vector<std::vector<std::vector<double>>> Polytopes::decomposePolytope (
    const std::vector<std::vector<double>>& vertices) const {
    size_t dimensions = vertices[0].size(); // Dimensionality (15 for 15D)
    size_t numPoints = vertices.size();

    // Flatten the vertices for Qhull
    std::vector<double> flatVertices;
    for (const auto& vertex : vertices) {
        flatVertices.insert(flatVertices.end(), vertex.begin(), vertex.end());
    }

    // Run Qhull to compute the convex hull
    orgQhull::Qhull qhull;
    qhull.runQhull("", dimensions, numPoints, flatVertices.data(), "Qt");

    // Compute the centroid of the polytope
    std::vector<double> centroid(dimensions, 0.0);
    for (const auto& vertex : vertices) {
        for (size_t i = 0; i < dimensions; ++i) {
            centroid[i] += vertex[i];
        }
    }
    for (size_t i = 0; i < dimensions; ++i) {
        centroid[i] /= numPoints;
    }

   
    std::vector<std::vector<std::vector<double>>> simplexes;
    for (const auto& facet : qhull.facetList()) {
        if (!facet.isGood()) continue; // Skip invalid facets

        // Extract facet vertices
        std::vector<std::vector<double>> simplex;
        for (const auto& vertex : facet.vertices()) {
            simplex.push_back(vertices[vertex.point().id()]);
        }

        // Add the centroid to form a 15D simplex
        simplex.push_back(centroid);

        // Validate simplex size
        if (simplex.size() != dimensions + 1) {
            std::cerr << "Warning: Invalid simplex size after adding centroid.\n";
            continue;
        }

        simplexes.push_back(simplex);
    }

    return simplexes;
}

// Calculate the volume of a simplex
double Polytopes::calculateSimplexVolume(const std::vector<std::vector<double>>& simplexVertices) const {
    size_t dimensions = simplexVertices[0].size(); 

    // Ensure the simplex has the correct number of vertices
    if (simplexVertices.size() != dimensions + 1) {
        std::cerr << "Error: Simplex does not have the correct number of vertices.\n";
        return 0.0;
    }

    // Construct the determinant matrix
    Eigen::MatrixXd matrix(dimensions, dimensions);
    for (size_t i = 1; i <= dimensions; ++i) {
        for (size_t j = 0; j < dimensions; ++j) {
            matrix(i - 1, j) = simplexVertices[i][j] - simplexVertices[0][j];
        }
    }

    // Compute the determinant
    double determinant = matrix.determinant();

    // Compute and return the volume
    return std::abs(determinant) / std::tgamma(dimensions + 1); // n! = Gamma(n+1)
}
