#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <vector>
#include <utility>
#include "../domain/integration_domain.h"


   // The Polytopes class represents a polytope (a geometric object with flat sides) in a multidimensional space.
// It extends IntegrationDomain, making it compatible with other domain-related functionalities.
class Polytopes : public IntegrationDomain {
 public:
    
    // Constructor: Initializes the polytope with its vertices.
    // Each vertex is represented as a vector of doubles.
   Polytopes(const std::vector<std::vector<double>>& vertices);
 
  // Determines if a given point is inside the polytope.
    // A point is represented as a vector of doubles.
    virtual bool contains(const std::vector<double>& point) const override;

   // Returns the number of dimensions of the polytope (e.g., 2D, 3D, etc.).
    size_t getDimensions() const;

    // Computes the bounded volume of the polytope.
    // This method will depend on the specific geometry of the polytope.
    virtual double getBoundedVolume() const override;

      // Decomposes the polytope into smaller simplices (e.g., triangles or tetrahedra) 
    // to facilitate volume computation or other operations.
    // Input: Vertices of the polytope.
    // Output: A collection of simplices, each represented as a set of vertices.
      std::vector<std::vector<std::vector<double>>> decomposePolytope(const std::vector<std::vector<double>> &vertices) const;

     // Calculates the volume of a single simplex (e.g., a triangle or tetrahedron) 
    // given its vertices.
    // Used as a building block for more complex volume calculations.
     double calculateSimplexVolume(const std::vector<std::vector<double>> &simplexVertices) const;

      // Provides the bounds of the polytope along each dimension as a collection of 
    // min-max pairs.
    virtual std::vector<std::pair<double, double>> getBounds() const override;

   
private:
      // Stores the vertices of the polytope.
    // Each vertex is a point in n-dimensional space.
     std::vector<std::vector<double>> vertices;          
};

#endif // POLYTOPES_H
