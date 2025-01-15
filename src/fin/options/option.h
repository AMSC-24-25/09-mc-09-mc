// option.h
#ifndef OPTION_H
#define OPTION_H

class Option {
public:
    enum class Type { Call, Put };

protected:
    Type optionType;
    double strike;
    double expiry;
    double spot;
    double volatility;
    double riskFreeRate;

public:
    Option(Type type, double strike, double expiry, double spot, double volatility, double riskFreeRate)
        : optionType(type), strike(strike), expiry(expiry), spot(spot), volatility(volatility), riskFreeRate(riskFreeRate) {}

    virtual ~Option() = default;
    virtual double payoff(double spotAtMaturity) const = 0;
    double getExpiry() const { return expiry; }
    double getSpot() const { return spot; }
    double getVolatility() const { return volatility; }
    double getRiskFreeRate() const { return riskFreeRate; }
};

#endif // OPTION_H
