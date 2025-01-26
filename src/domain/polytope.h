#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <vector>
#include <utility>
#include "../domain/integration_domain.h"

class Polytopes : public IntegrationDomain {
public:
    /**
     * @brief Constructs a polytope defined by \( Ax \leq b \).
     * 
     * @param A The matrix of coefficients where each row corresponds to a constraint.
     * @param b The vector of bounds corresponding to each constraint.
     */
    Polytopes(const std::vector<std::vector<double>>& A, const std::vector<double>& b);

    /**
     * @brief Checks if a point lies within the polytope.
     * 
     * A point is considered inside if it satisfies \( Ax \leq b \).
     *
     * @param point The point to test.
     * @return true if the point is inside the polytope.
     * @return false otherwise.
     */
    virtual bool contains(const std::vector<double>& point) const override;

    /**
     * @brief Returns the dimensions of the polytope.
     *
     * @return size_t The number of dimensions.
     */
    size_t getDimensions() const;

    /**
     * @brief Calculates the volume of the polytope.
     * 
     * This is an optional method for debugging purposes or for specific use cases
     * where an approximate volume is needed.
     * 
     * @return double Approximate volume (default: -1 for no calculation).
     */
    virtual double getBoundedVolume() const override;

    /**
     * @brief Returns an axis-aligned bounding box for the polytope.
     * 
     * @return std::vector<std::pair<double, double>> The bounds for each dimension.
     */
    virtual std::vector<std::pair<double, double>> getBounds() const override;



private:
    std::vector<std::vector<double>> A;  // Coefficient matrix for the inequalities
    std::vector<double> b;               // Right-hand side values of the inequalities
};

#endif // POLYTOPES_H
