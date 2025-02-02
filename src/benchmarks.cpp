#include "benchmarks.h"

const std::vector<size_t> numPointsValues = {10'000, 100'000, 1'000'000, 10'000'000, 100'000'000, 1'000'000'000};
unsigned int numThreads;

//Function to export results in txt
void exportIntegrationResults(const std::string &fileName,
                             const std::vector<MCResultRow> &rows)
{
    std::ofstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Impossibile aprire il file " << fileName << " per scrittura.\n";
        return;
    }

    // File Header here
    file << "NumPoints\tGridDim\tTimeStd(ms)\tTimeStrat(ms)\tStdResult\tStratResult\n";

    // write the rows
    for (const auto &row : rows) {
        file << row.numPoints    << "\t"
             << row.gridDim      << "\t"
             << row.timeStd      << "\t"
             << row.timeStrat    << "\t"
             << row.stdResult    << "\t"
             << row.stratResult  << "\n";
    }

    file.close();
}


void circleIntegration() {
    //Function to integrate on the circle: x^2 + y^2
    auto f = [](const std::vector<double> &x) {
        return x[0] * x[0] + x[1] * x[1];
    };

    // Circle with radius 1 centered at (0, 0)
    Hypersphere sphere(2, 1.0);

    std::vector<MCResultRow> circleData; //Vector to store rows

    // Monte Carlo integrator initialization
    MonteCarloIntegrator mcIntegrator(sphere);

    // Print header
    std::cout << "Integrating f(x,y) = x^2 + y^2 over the unit circle (radius = 1):\n";
    std::cout << "Expected result (π/2): " << std::fixed << std::setprecision(6) << (M_PI / 2) << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Time Strat (ms)"
              << std::setw(18) << "Std Result"
              << std::setw(20) << "Strat Result" << "\n";
    std::cout << std::string(98, '-') << "\n";

    //Different sizes for the layered method
    std::vector<int32_t> strataPerDimValues = {5, 10, 20, 50};
    //std::vector<int32_t> strataPerDimValues = {10, 100, 1000, 5000};        //trying something different

    for (size_t numPoints : numPointsValues) {
        // Standard method
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);

        //Saving results in the struct
         {
        MCResultRow row;
        row.numPoints    = numPoints;
        row.gridDim      = "N/A";

        long long timeStdMs = durationStandard.count();
        row.timeStd   = std::to_string(timeStdMs);
        row.timeStrat    = "-";
        // Convert results to string
        std::ostringstream ossRes;
        ossRes << std::fixed << std::setprecision(6) << resultStandard;
        row.stdResult    = ossRes.str();
        row.stratResult  = "-";

        circleData.push_back(row);
    }
        // Printing the standard method line
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"  //No grid for the standard method
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << "-"  // No time for the layered method
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard
                  << std::setw(20) << "-" << "\n";

        // Stratified method
        for (int32_t strataPerDim : strataPerDimValues) {
            auto startStrat = std::chrono::high_resolution_clock::now();
            double resultStratified = mcIntegrator.integrateStratified(f, numPoints, numThreads, strataPerDim);
            auto endStrat = std::chrono::high_resolution_clock::now();
            auto durationStratified = std::chrono::duration_cast<std::chrono::milliseconds>(endStrat - startStrat);

        //Saving results in the struct
         {
            MCResultRow row;
            row.numPoints   = numPoints;
            row.gridDim     = std::to_string(strataPerDim);
            row.timeStd     = "-";
            long long timeStratMs = durationStratified.count();
            row.timeStrat = std::to_string(timeStratMs);
            row.stdResult   = "-";

            std::ostringstream ossRes;
            ossRes << std::fixed << std::setprecision(6) << resultStratified;
            row.stratResult = ossRes.str();

            circleData.push_back(row);
        }

            // Printing of the layered method line
            std::cout << std::setw(12) << numPoints
                      << std::setw(12) << strataPerDim
                      << std::setw(18) << "-"  // No time for the standard method
                      << std::setw(18) << durationStratified.count()
                      << std::setw(18) << "-"
                      << std::setw(20) << resultStratified << "\n";
        }
    }
    std::cout << "\nSaved txt file: resultsCircle.txt\n";
    exportIntegrationResults("resultsCircle.txt", circleData);
}


void triangleIntegration() {
    // Integratiion over the (1,1,1) triangle domain with f(x,y)=1
    // to find the area of this triangle.
    // The vertices of an equilateral triangle of side length 1 are as follows:
    // A=(0,0), B=(1,0), C=(0.5, sqrt(3)/2)
    std::vector<std::pair<double, double>> triangleVertices = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.5, std::sqrt(3.0) / 2.0}
    };
    Polygon2D triangle(triangleVertices);

    std::vector<MCResultRow> triangleData; //Vector to store rows

    // Monte Carlo integrator for the triangle
    MonteCarloIntegrator mcIntegrator(triangle);

    // The function to integrate is f(x,y)=1
    auto f = [](const auto &) {
        return 1.0;
    };

    // The expected result is the area of the triangle = sqrt(3)/4
    double expectedTriangleArea = std::sqrt(3.0) / 4.0;

    // Print header
    std::cout << "\nIntegrating f(x,y)=1 over the equilateral triangle (1,1,1):\n";
    std::cout << "Expected area: " << std::fixed << std::setprecision(6) << expectedTriangleArea << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Time Strat (ms)"
              << std::setw(18) << "Std Result"
              << std::setw(20) << "Strat Result" << "\n";
    std::cout << std::string(98, '-') << "\n";

    // Different sizes for the layered method
    std::vector<int32_t> strataPerDimValues = {5, 10, 20, 50};
    //std::vector<int32_t> strataPerDimValues = {10, 100, 1000, 5000};     //trying something different

    for (size_t numPoints : numPointsValues) {
        // Standard method
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);


        //Saving results in the struct
         {
        MCResultRow row;
        row.numPoints    = numPoints;
        row.gridDim      = "N/A";

        long long timeStdMs = durationStandard.count();
        row.timeStd   = std::to_string(timeStdMs);
        row.timeStrat    = "-";
        // Convert results to string
        std::ostringstream ossRes;
        ossRes << std::fixed << std::setprecision(6) << resultStandard;
        row.stdResult    = ossRes.str();
        row.stratResult  = "-";

        triangleData.push_back(row);
    }

        // Printing the line for the standard method
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"                     // No grid for the standard method
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << "-"                       // No time for the layered method
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard
                  << std::setw(20) << "-" << "\n";

        // Stratified Method
        for (int32_t strataPerDim : strataPerDimValues) {
            auto startStrat = std::chrono::high_resolution_clock::now();
            double resultStratified = mcIntegrator.integrateStratified(f, numPoints, numThreads, strataPerDim);
            auto endStrat = std::chrono::high_resolution_clock::now();
            auto durationStratified = std::chrono::duration_cast<std::chrono::milliseconds>(endStrat - startStrat);


            //Save results in the struct
            {
            MCResultRow row;
            row.numPoints   = numPoints;
            row.gridDim     = std::to_string(strataPerDim);
            row.timeStd     = "-";
            long long timeStratMs = durationStratified.count();
            row.timeStrat = std::to_string(timeStratMs);
            row.stdResult   = "-";

            std::ostringstream ossRes;
            ossRes << std::fixed << std::setprecision(6) << resultStratified;
            row.stratResult = ossRes.str();

            triangleData.push_back(row);
        }

            // Printing the row for the layered method
            std::cout << std::setw(12) << numPoints
                      << std::setw(12) << strataPerDim
                      << std::setw(18) << "-"             // No time for the standard method
                      << std::setw(18) << durationStratified.count()
                      << std::setw(18) << "-"
                      << std::setw(20) << std::fixed << std::setprecision(6) << resultStratified << "\n";
        }
    }
    std::cout << "\nSaved txt file: resultsTriangle.txt\n";
    exportIntegrationResults("resultsTriangle.txt", triangleData);
}

void fiveDimIntegration() {

    // Define f in 5 dim
    //    f(x0,x1,x2,x3,x4) = x0^2 + x1^2 + x2^2 + x3^2 + x4^2
    auto f = [](const std::vector<double> &x) {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4];
    };

    // Hyperspher r=1
    Hypersphere sphere(5, 1.0);

    std::vector<MCResultRow> Function5dData; //Vector to store rows

    MonteCarloIntegrator mcIntegrator(sphere);

    // Print
    std::cout << "\nIntegrating f(x,y,z,u,w) = x^2 + y^2 + z^2 + u^2 + w^2 "
              << "over the 5D unit hypersphere (radius = 1):\n";
    std::cout << "Expected result (8π^2/21) ~ 3.75985 " << "\n\n";
    

    // Print header
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Time Strat (ms)"
              << std::setw(18) << "Std Result"
              << std::setw(20) << "Strat Result" << "\n";
    std::cout << std::string(98, '-') << "\n";

    // Choose stratadim (strataPerDim),
    std::vector<int32_t> strataPerDimValues = {5, 10, 20, 50};

    
    for (size_t numPoints : numPointsValues) {
        // ------ Standard method ------
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);


        //Saving results in the struct
         {
        MCResultRow row;
        row.numPoints    = numPoints;
        row.gridDim      = "N/A";
        
        long long timeStdMs = durationStandard.count();
        row.timeStd   = std::to_string(timeStdMs);
        row.timeStrat    = "-";   
        // Convert results to string
        std::ostringstream ossRes;
        ossRes << std::fixed << std::setprecision(6) << resultStandard;
        row.stdResult    = ossRes.str();
        row.stratResult  = "-";
        
        Function5dData.push_back(row);
    }

        // Print results
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"  // No grid for std
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << "-"  // No time for strat
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard
                  << std::setw(20) << "-" << "\n";
       
        // ------ Stratified ------
        for (int32_t strataPerDim : strataPerDimValues) {
            auto startStrat = std::chrono::high_resolution_clock::now();
            double resultStratified = mcIntegrator.integrateStratified(f, numPoints, numThreads, strataPerDim);
            auto endStrat = std::chrono::high_resolution_clock::now();
            auto durationStratified = std::chrono::duration_cast<std::chrono::milliseconds>(endStrat - startStrat);

        //Saving results in the struct
        {
            MCResultRow row;
            row.numPoints   = numPoints;
            row.gridDim     = std::to_string(strataPerDim);
            row.timeStd     = "-";
            long long timeStratMs = durationStratified.count();
            row.timeStrat = std::to_string(timeStratMs);
            row.stdResult   = "-";
            
            std::ostringstream ossRes;
            ossRes << std::fixed << std::setprecision(6) << resultStratified;
            row.stratResult = ossRes.str();

            Function5dData.push_back(row);
        }


            // Print strat results
            std::cout << std::setw(12) << numPoints
                      << std::setw(12) << strataPerDim
                      << std::setw(18) << "-"  // No std time
                      << std::setw(18) << durationStratified.count()
                      << std::setw(18) << "-" 
                      << std::setw(20) << std::fixed << std::setprecision(6) << resultStratified
                      << "\n";

                     
        }
    }

    
    exportIntegrationResults("resultsFunction5D.txt", Function5dData);
    std::cout << "\nSaved txt file: resultsFunction5D.txt\n";

}

void twelveDimIntegration() {
    constexpr int dimensions = 12;

    // Define f in 12 dim
    //    f(x) = sum_{i=0}^11 x_i^2
    auto f = [](const std::vector<double> &x) {
        double res = 0;
        for (size_t i = 0; i < dimensions; i++) {
            res += x[i] * x[i];
        }
        return res;
    };

    // Hypersphere r=1
    Hypersphere sphere(dimensions, 1.0);

    std::vector<MCResultRow> Function12dData; // Vector to store rows

    MonteCarloIntegrator mcIntegrator(sphere);

    // Print
    std::cout << "\nIntegrating f(x) = sum_{i=0}^11 x_i^2 "
              << "over the area inside the 12D unit hypersphere (radius = 1):\n";
    std::cout << "Expected result (π^6/840) ~ 1.14451 " << "\n\n";


    // Print header
    std::cout << std::setw(12) << "NumPoints"
              << std::setw(12) << "GridDim"
              << std::setw(18) << "Time Std (ms)"
              << std::setw(18) << "Std Result" << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t numPoints : numPointsValues) {
        // ------ Standard method ------
        auto startStd = std::chrono::high_resolution_clock::now();
        double resultStandard = mcIntegrator.integrate(f, numPoints, numThreads);
        auto endStd = std::chrono::high_resolution_clock::now();
        auto durationStandard = std::chrono::duration_cast<std::chrono::milliseconds>(endStd - startStd);


        //Saving results in the struct
        {
            MCResultRow row;
            row.numPoints = numPoints;
            row.gridDim = "N/A";

            long long timeStdMs = durationStandard.count();
            row.timeStd = std::to_string(timeStdMs);
            row.timeStrat = "-";
            // Convert results to string
            std::ostringstream ossRes;
            ossRes << std::fixed << std::setprecision(6) << resultStandard;
            row.stdResult = ossRes.str();
            row.stratResult = "-";

            Function12dData.push_back(row);
        }

        // Print results
        std::cout << std::setw(12) << numPoints
                  << std::setw(12) << "N/A"  // No grid for std
                  << std::setw(18) << durationStandard.count()
                  << std::setw(18) << std::fixed << std::setprecision(6) << resultStandard << "\n";
    }

    exportIntegrationResults("resultsFunction12D.txt", Function12dData);
    std::cout << "\nSaved txt file: resultsFunction12D.txt\n";
}


void MHintegration() {
    // Suppose we have a standard 2D normal distribution
    auto p = [](const std::vector<double> &x) -> double {
        return 1.0 / (2 * M_PI) * exp(-(x[0] * x[0] + x[1] * x[1]) / 2);
    };

    // Let's calculate E[x_0^2]
    // It should be equal to 1.0 for standard normal distribution
    auto f = [](const std::vector<double> &x) -> double {
        return x[0] * x[0];
    };

    // Integrate the square with the most probability density
    Polygon2D domain({
        {-10.0, -10.0},
        {-10.0, 10.0},
        {10.0, 10.0},
        {10.0, -10.0}
    });

    // MH integrator for the square (set stddev for the Gaussian noise to 0.5)
    MetropolisHastingsIntegrator mh(domain, 0.5);
    std::vector<double> initialPoint = {0.0, 0.0};

    std::cout << "\nCalculating E[x_0^2] of a standard 2D normal distribution:\n";
    std::cout << "Expected result: " << std::fixed << std::setprecision(4) << 1.0 << "\n\n";
    std::cout << std::setw(12) << "NumPoints"
            << std::setw(12) << "Time (ms)"
            << std::setw(13) << "Result"
            << std::setw(18) << "Acceptance rate" << "\n";
    std::cout << std::string(83, '-') << "\n";

    for (size_t numPoints: numPointsValues) {
        // MH is more computationally-heavy
        size_t points = numPoints / 10;

        // Metropolis-Hastings integration on the square
        auto start = std::chrono::high_resolution_clock::now();
        auto [result, acceptanceRate] = mh.integrateParallel(f, p, points, initialPoint, numThreads);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << std::setw(12) << points
                << std::setw(11) << duration.count()
                << std::setw(14) << std::fixed << std::setprecision(6) << result
                << std::setw(18) << std::fixed << std::setprecision(2) << acceptanceRate * 100 << "%\n";
    }
}

void tetrahedron_integration(){

    std::vector<std::vector<double>> tetrahedronVertices = {
        {0.0, 0.0, 0.0},  // Vertex 1
        {1.0, 0.0, 0.0},  // Vertex 2
        {0.0, 1.0, 0.0},  // Vertex 3
        {0.0, 0.0, 1.0}   // Vertex 4
    };

    Polytopes tetrahedron(tetrahedronVertices);
    MonteCarloIntegrator mcIntegrator(tetrahedron);

    auto f = [](const std::vector<double>& point) -> double {
        // f(x,y,z) = x + y
        return point[0] + point[1];
    };

      std::cout << "\nIntegrating f(x,y) = x + y "
              << "over the area inside the 3D unit hypercube :\n";
    std::cout << "Expected result 0.0833333 " << "\n\n";

    std::vector<MCResultRow> tetraData;
    for (size_t numPoints : numPointsValues) {
        auto start = std::chrono::high_resolution_clock::now();
        double result = mcIntegrator.integrate(f, numPoints, numThreads);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        MCResultRow row;
        row.numPoints = numPoints;
        row.gridDim   = "N/A";
        row.timeStd   = std::to_string(duration.count());
        row.timeStrat = "-";
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << result;
        row.stdResult    = oss.str();
        row.stratResult  = "-";
        tetraData.push_back(row);

    
        std::cout << "Tetrahedron: NumPoints=" << numPoints 
                  << " Time=" << duration.count() << " ms "
                  << " Result=" << result << "\n";
    }
    exportIntegrationResults("resultsTetrahedron.txt", tetraData);
    std::cout << "\nSaved txt file: resultsTetrahedron.txt\n";
}


void simplex5DIntegration() {
    std::vector<std::vector<double>> simplex5D = {
        {0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0}
    };

    Polytopes simplex(simplex5D);
    MonteCarloIntegrator mcIntegrator(simplex);

    auto f = [](const std::vector<double>& point) -> double {
        // f(x0,x1,x2,x3,x4)= x0+x1+x2
        return point[0] + point[1] + point[2];
    };

     std::cout << "\nIntegrating f = x + y + z "
              << "over the area inside the 5D unit hypercube :\n";
    std::cout << "Expected result 0.00417 " << "\n\n";

    std::vector<MCResultRow> simplexData;
    for (size_t numPoints : numPointsValues) {
        auto start = std::chrono::high_resolution_clock::now();
        double result = mcIntegrator.integrate(f, numPoints, numThreads);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        MCResultRow row;
        row.numPoints = numPoints;
        row.gridDim   = "N/A";
        row.timeStd   = std::to_string(duration.count());
        row.timeStrat = "-";
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << result;
        row.stdResult    = oss.str();
        row.stratResult  = "-";
        simplexData.push_back(row);

        std::cout << "5D Simplex: NumPoints=" << numPoints 
                  << " Time=" << duration.count() << " ms "
                  << " Result=" << result << "\n";
    }
    exportIntegrationResults("resultsSimplex5D.txt", simplexData);
    std::cout << "\nSaved txt file: resultsSimplex5D.txt\n";
}




void benchmarks() {
    numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) {
        // fallback
        numThreads = 16;
    }
    std::cout << "Using " << numThreads << " threads.\n";

   circleIntegration();
   triangleIntegration();
   fiveDimIntegration();
   twelveDimIntegration();
   MHintegration();
    tetrahedron_integration();
     simplex5DIntegration();

}
