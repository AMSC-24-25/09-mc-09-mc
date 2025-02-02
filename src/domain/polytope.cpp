#include "polytope.h"
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullVertexSet.h> 
#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>
#include <numeric>

// Constructor
Polytopes::Polytopes(const std::vector<std::vector<double>>& vertices)
    : vertices(vertices) {
    if (vertices.empty()) {
        throw std::invalid_argument("The vertices cannot be empty.");
    }
}

// This version works for any n-dimensional simplex (n+1 vertices in R^n).
bool Polytopes:: contains(const std::vector<double>& point) const {
        size_t n_points = vertices.size();
        size_t n_dim = point.size();

        if (n_dim != vertices[0].size())
            throw std::invalid_argument("Point dimension does not match polytope dimension.");

        // For an n-dimensional simplex, we require exactly n+1 vertices.
        if (n_points != n_dim + 1)
            throw std::invalid_argument("The polytope is not an n-dimensional simplex (expected n+1 vertices).");

        // Use the first vertex as the reference: v0.
        Eigen::VectorXd v0(n_dim);
        for (size_t j = 0; j < n_dim; ++j)
            v0(j) = vertices[0][j];

        // Build the n x n matrix M whose columns are: v_i - v0 for i = 1 ... n.
        Eigen::MatrixXd M(n_dim, n_dim);
        for (size_t i = 0; i < n_dim; ++i) {
            for (size_t j = 0; j < n_dim; ++j) {
                M(j, i) = vertices[i + 1][j] - v0(j);
            }
        }

        // Compute the right-hand side: p - v0, where p is the given point.
        Eigen::VectorXd rhs(n_dim);
        for (size_t j = 0; j < n_dim; ++j)
            rhs(j) = point[j] - v0(j);

        // Solve for the barycentric coordinates corresponding to vertices 1...n.
        // (These are the λ_i for i=1...n.)
        Eigen::VectorXd lambda = M.colPivHouseholderQr().solve(rhs);

        // Compute lambda0 = 1 - (λ_1 + λ_2 + ... + λ_n)
        double lambda0 = 1.0 - lambda.sum();

        // Allow a small numerical tolerance.
        double tol = 1e-6;

        // For the point to be inside, each barycentric coordinate must be >= -tol.
        if (lambda0 < -tol)
            return false;
        for (int i = 0; i < lambda.size(); ++i) {
            if (lambda(i) < -tol)
                return false;
        }
        return true;
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
    auto bounds = getBounds();
        double volume = 1.0;
        for (const auto& b : bounds)
            volume *= (b.second - b.first);
        return volume;
}

// Function to decompose a nD polytope into simplexes
std::vector<std::vector<std::vector<double>>> Polytopes::decomposePolytope (
    const std::vector<std::vector<double>>& vertices) const {
    size_t dimensions = vertices[0].size(); // Dimensionality 
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

        // Add the centroid to form a nD simplex
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
