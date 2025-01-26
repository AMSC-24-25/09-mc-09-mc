#ifndef POLYTOPES_H
#define POLYTOPES_H

#include <vector>
#include <utility>
#include "../domain/integration_domain.h"

class Polytopes : public IntegrationDomain {
public:
   
   // Polytopes(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
   Polytopes(const std::vector<std::vector<double>>& vertices);

    virtual bool contains(const std::vector<double>& point) const override;
    size_t getDimensions() const;
    virtual double getBoundedVolume() const override;
      std::vector<std::vector<std::vector<double>>> decomposePolytope(const std::vector<std::vector<double>> &vertices) const;

     double calculateSimplexVolume(const std::vector<std::vector<double>> &simplexVertices) const;

    virtual std::vector<std::pair<double, double>> getBounds() const override;

private:
    private:
    std::vector<std::vector<double>> vertices;          
};

#endif // POLYTOPES_H
