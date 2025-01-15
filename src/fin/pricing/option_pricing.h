#ifndef OPTION_PRICING_H
#define OPTION_PRICING_H

#include "../options/option.h"
#include "../../integrators/abstract_integrator.h"

class OptionPricing {
public:
    static double priceOption(const Option &option, size_t numSimulations, AbstractIntegrator &integrator);
};

#endif // OPTION_PRICING_H
