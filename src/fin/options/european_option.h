// european_option.h
#ifndef EUROPEAN_OPTION_H
#define EUROPEAN_OPTION_H

#include "option.h"

class EuropeanOption : public Option {
public:
    EuropeanOption(Type type, double strike, double expiry, double spot, double volatility, double riskFreeRate)
        : Option(type, strike, expiry, spot, volatility, riskFreeRate) {}

    double payoff(double spotAtMaturity) const override {
        if (optionType == Type::Call) {
            return std::max(spotAtMaturity - strike, 0.0);
        } else {
            return std::max(strike - spotAtMaturity, 0.0);
        }
    }
};

#endif // EUROPEAN_OPTION_H
