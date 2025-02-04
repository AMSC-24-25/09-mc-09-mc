#include "room_assignment.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

RoomAssignment::RoomAssignment(int num_students, const std::vector<std::vector<double>>& dislike_matrix)
    : n(num_students), dislikes(dislike_matrix), gen(rd()) {
    if (n % 2 != 0) {
        throw std::invalid_argument("Number of students must be even");
    }
    assignments.resize(n);

    // Initialize random assignments
    for (int i = 0; i < n; i++) {
        assignments[i] = i / 2;  // Initial assignment: pairs of students
    }

    // Shuffle the assignments randomly
    std::shuffle(assignments.begin(), assignments.end(), gen);
}

double RoomAssignment::calculateTotalDislikes() {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (assignments[i] == assignments[j] && i != j) {
                sum += dislikes[i][j];
            }
        }
    }
    return sum;
}

std::pair<int, int> RoomAssignment::selectRandomStudents() {//@note this method should be const
    std::uniform_int_distribution<> dis(0, n - 1);
    int s1 = dis(gen);
    int s2;
    do {
        s2 = dis(gen);
    } while (assignments[s1] == assignments[s2]);
    return {s1, s2};
}

void RoomAssignment::recordMetrics(int iteration, double temperature, double cost, double acceptance_rate) {
    metrics.push_back({
        iteration,
        temperature,
        cost,
        acceptance_rate
    });
}

void RoomAssignment::solve() {
    double temperature = 1.0;//@note better put this stuff in a struct, maybe with a method that anable to change the defaults
    int iterations = 0;
    int max_iterations = 1000;
    int no_changes = 0;
    int accepted_moves = 0;

    double current_cost = calculateTotalDislikes();
    int metric_interval = max_iterations / 5;

    while (no_changes < max_iterations && temperature > 0.01) {
        auto [student1, student2] = selectRandomStudents();

        // Store original assignments
        int room1 = assignments[student1];
        int room2 = assignments[student2];

        // Swap rooms
        assignments[student1] = room2;
        assignments[student2] = room1;

        double new_cost = calculateTotalDislikes();
        double delta = new_cost - current_cost;

        // Accept or reject new solution
        std::uniform_real_distribution<> dis(0.0, 1.0);
        if (delta < 0 || dis(gen) < exp(-delta / temperature)) {
            current_cost = new_cost;
            no_changes = 0;
            accepted_moves++;
        } else {
            // Revert swap
            assignments[student1] = room1;
            assignments[student2] = room2;
            no_changes++;
        }

        if (iterations % metric_interval == 0) {
            double acceptance_rate = static_cast<double>(accepted_moves) / (iterations + 1);
            recordMetrics(iterations, temperature, current_cost, acceptance_rate);
        }

        temperature *= 0.999;//@note better not have contants hardwired. Treat it as a parameter
        iterations++;
    }

    // Record final state
    double final_acceptance_rate = static_cast<double>(accepted_moves) / iterations;
    recordMetrics(iterations, temperature, current_cost, final_acceptance_rate);
}

void RoomAssignment::printMetricsTable() const {
    // Print header
    std::cout << std::setw(12) << "Iteration"
              << std::setw(15) << "Temperature"
              << std::setw(15) << "Cost"
              << std::setw(15) << "Accept Rate"
              << "\n";
    std::cout << std::string(57, '-') << "\n";

    // Print metrics
    for (const auto& metric : metrics) {
        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(12) << metric.iteration
                  << std::setw(15) << metric.temperature
                  << std::setw(15) << metric.cost
                  << std::setw(15) << metric.acceptance_rate
                  << "\n";
    }
}

void RoomAssignment::printAssignments() {
    std::cout << "Room assignments:\n";
    for (int room = 0; room < n/2; room++) {
        std::cout << "Room " << room << ": ";
        for (int student = 0; student < n; student++) {
            if (assignments[student] == room) {
                std::cout << student << " ";
            }
        }
        std::cout << "\n";
    }
}

//@note I would have separated the tools of your library from the examples
void room_assignment_example() {
    // Example with 12 students (6 rooms)
    int n = 12;

    // A complex dislike matrix for 12 students
    // Values range from 0.0 (no dislike) to 1.0 (maximum dislike)
    std::vector<std::vector<double>> dislikes = {
        //  0    1    2    3    4    5    6    7    8    9    10   11
        {0.0, 0.8, 0.4, 0.3, 0.2, 0.1, 0.5, 0.6, 0.7, 0.3, 0.2, 0.4}, // Student 0
        {0.8, 0.0, 0.3, 0.6, 0.4, 0.2, 0.3, 0.4, 0.8, 0.5, 0.1, 0.3}, // Student 1
        {0.4, 0.3, 0.0, 0.5, 0.9, 0.4, 0.2, 0.3, 0.4, 0.6, 0.7, 0.5}, // Student 2
        {0.3, 0.6, 0.5, 0.0, 0.3, 0.8, 0.4, 0.5, 0.2, 0.4, 0.6, 0.8}, // Student 3
        {0.2, 0.4, 0.9, 0.3, 0.0, 0.5, 0.6, 0.7, 0.3, 0.2, 0.4, 0.6}, // Student 4
        {0.1, 0.2, 0.4, 0.8, 0.5, 0.0, 0.7, 0.8, 0.4, 0.3, 0.5, 0.2}, // Student 5
        {0.5, 0.3, 0.2, 0.4, 0.6, 0.7, 0.0, 0.9, 0.5, 0.4, 0.3, 0.1}, // Student 6
        {0.6, 0.4, 0.3, 0.5, 0.7, 0.8, 0.9, 0.0, 0.6, 0.5, 0.2, 0.3}, // Student 7
        {0.7, 0.8, 0.4, 0.2, 0.3, 0.4, 0.5, 0.6, 0.0, 0.7, 0.4, 0.5}, // Student 8
        {0.3, 0.5, 0.6, 0.4, 0.2, 0.3, 0.4, 0.5, 0.7, 0.0, 0.8, 0.6}, // Student 9
        {0.2, 0.1, 0.7, 0.6, 0.4, 0.5, 0.3, 0.2, 0.4, 0.8, 0.0, 0.7}, // Student 10
        {0.4, 0.3, 0.5, 0.8, 0.6, 0.2, 0.1, 0.3, 0.5, 0.6, 0.7, 0.0}  // Student 11
    };

    try {//@note nice using exceptions
        // Create solver with 12 students
        RoomAssignment solver(n, dislikes);

        // Print initial state
        std::cout << "Initial arrangement:\n";
        solver.printAssignments();
        std::cout << "Initial cost: " << solver.calculateTotalDislikes() << "\n\n";

        // Run the simulated annealing algorithm
        std::cout << "Running simulated annealing optimization...\n";
        solver.solve();

        std::cout << "\nPerformance Metrics:\n";
        solver.printMetricsTable();

        // Print final state
        std::cout << "\nFinal arrangement:\n";
        solver.printAssignments();

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
}
