// benchmark.h

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "integrators/metropolis_hastings_integrator.h"
#include "integrators/monte_carlo_integrator.h"
#include "domain/hypersphere.h"
#include "domain/polygon2d.h"
#include "domain/polytope.h"

#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <cmath>

// Struct to store results in a txt file
struct MCResultRow {
    size_t numPoints;         
    std::string gridDim;      
    std::string timeStd;      
    std::string timeStrat;    
    std::string stdResult;    
    std::string stratResult;  
};

// Function prototypes
void exportIntegrationResults(const std::string &fileName, const std::vector<MCResultRow> &rows);
void circleIntegration();
void triangleIntegration();
void fiveDimIntegration();
void twelveDimIntegration();
void MHintegration();
void benchmarks();

#endif // BENCHMARK_H