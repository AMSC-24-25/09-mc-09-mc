#ifndef DUMMY_DOMAIN_H
#define DUMMY_DOMAIN_H

#include "integration_domain.h"
#include <vector>
#include <utility>

// Dummy domain class that implements IntegrationDomain interface
class DummyDomain : public IntegrationDomain {
public:
    DummyDomain() = default;
    // Returns empty bounds (or any dummy bounds)
    std::vector<std::pair<double, double>> getBounds() const override {
        return {};
    }

    // Returns a dummy volume (could be zero or one)
    double getBoundedVolume() const override {
        return 1.0;
    }

    // Always returns true, as the dummy domain contains all points
    bool contains(const std::vector<double> &point) const override {
        return true;
    }
};

#endif // DUMMY_DOMAIN_H
