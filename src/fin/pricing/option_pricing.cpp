#include "option_pricing.h"

double OptionPricing::priceOption(const Option &option, size_t numSimulations, AbstractIntegrator &integrator) {
    double sumPayoffs = 0.0;
    double dt = option.getExpiry();
    double spot = option.getSpot();
    double volatility = option.getVolatility();
    double riskFreeRate = option.getRiskFreeRate();
    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> dist(0.0, 1.0);

    for (size_t i = 0; i < numSimulations; ++i) {
        double gauss_bm = dist(gen);
        double spotAtMaturity = spot * exp((riskFreeRate - 0.5 * volatility * volatility) * dt + volatility * sqrt(dt) * gauss_bm);
        double payoff = option.payoff(spotAtMaturity);
        sumPayoffs += payoff;
    }

    return (sumPayoffs / static_cast<double>(numSimulations)) * exp(-riskFreeRate * dt); // Sconto del payoff atteso
}
