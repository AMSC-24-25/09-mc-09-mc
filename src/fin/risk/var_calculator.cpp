#include "var_calculator.h"
#include <algorithm>

double VarCalculator::calculateVaR(const std::vector<double> &portfolioValues, double confidenceLevel) {
    size_t index = static_cast<size_t>((1.0 - confidenceLevel) * portfolioValues.size());
    std::vector<double> sortedValues = portfolioValues;
    std::sort(sortedValues.begin(), sortedValues.end());
    return sortedValues[index];
}
