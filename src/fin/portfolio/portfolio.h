#ifndef PORTFOLIO_H
#define PORTFOLIO_H

#include <vector>
#include <memory>
#include "../options/option.h"
#include "../pricing/option_pricing.h"
#include "../../integrators/abstract_integrator.h"

class Portfolio {
private:
    std::vector<std::shared_ptr<Option>> options;

public:
    void addOption(const std::shared_ptr<Option>& option);
    double calculatePortfolioValue(size_t numSimulations, AbstractIntegrator &integrator) const;
};

#endif // PORTFOLIO_H
