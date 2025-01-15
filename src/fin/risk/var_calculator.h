// var_calculator.h
#ifndef VAR_CALCULATOR_H
#define VAR_CALCULATOR_H

#include <vector>
#include <algorithm>

class VarCalculator {
public:
    static double calculateVaR(const std::vector<double> &portfolioValues, double confidenceLevel);
};

#endif // VAR_CALCULATOR_H
