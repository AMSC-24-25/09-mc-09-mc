#include "portfolio.h"

void Portfolio::addOption(const std::shared_ptr<Option>& option) {
    options.push_back(option);
}

double Portfolio::calculatePortfolioValue(size_t numSimulations, AbstractIntegrator &integrator) const {
    double portfolioValue = 0.0;
    for (const auto& option : options) {
        portfolioValue += OptionPricing::priceOption(*option, numSimulations, integrator);
    }
    return portfolioValue;
}
