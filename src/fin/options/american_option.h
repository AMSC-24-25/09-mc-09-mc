/* #ifndef AMERICAN_OPTION_H
#define AMERICAN_OPTION_H

#include "option.h"
#include <vector> // for std::vector
#include <algorithm> // for std::max
#include <cmath> // for std::exp, std::sqrt

class AmericanOption : public Option {
public:
    AmericanOption(Type type, double strike, double expiry, double spot, double volatility, double riskFreeRate)
        : Option(type, strike, expiry, spot, volatility, riskFreeRate) {}

    double payoff(double spotAtMaturity) const override {
        if (optionType == Type::Call) {
            return std::max(spotAtMaturity - strike, 0.0);
        } else {
            return std::max(strike - spotAtMaturity, 0.0);
        }
    }

    double price(int steps) const {
        std::vector<std::vector<double>> priceTree(steps + 1, std::vector<double>(steps + 1));
        std::vector<std::vector<double>> valueTree(steps + 1, std::vector<double>(steps + 1));
        
        double deltaT = expiry / steps;
        double discount = std::exp(-riskFreeRate * deltaT);
        double u = std::exp(volatility * std::sqrt(deltaT));
        double d = 1.0 / u;
        double pu = (std::exp(riskFreeRate * deltaT) - d) / (u - d);
        double pd = 1.0 - pu;

        // Initialize the price tree
        for (int i = 0; i <= steps; ++i) {
            for (int j = 0; j <= i; ++j) {
                priceTree[i][j] = spot * std::pow(u, j) * std::pow(d, i - j);
            }
        }

        // Initialize the option value at maturity
        for (int j = 0; j <= steps; ++j) {
            valueTree[steps][j] = payoff(priceTree[steps][j]);
        }

        // Backward induction
        for (int i = steps - 1; i >= 0; --i) {
            for (int j = 0; j <= i; ++j) {
                double exerciseValue = payoff(priceTree[i][j]);
                double continuationValue = discount * (pu * valueTree[i + 1][j + 1] + pd * valueTree[i + 1][j]);
                valueTree[i][j] = std::max(exerciseValue, continuationValue);
            }
        }

        return valueTree[0][0];
    }
};

#endif // AMERICAN_OPTION_H
*/