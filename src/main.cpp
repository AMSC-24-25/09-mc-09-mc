#include "benchmarks.h"
#include "ising_model.h"
#include "room_assignment.h"

#include <iostream>

void showMenu() {
    std::cout << "\nWhat would you like to run?\n";
    std::cout << "1. Benchmarks for each integration and domain\n";
    std::cout << "2. Specific heat per particle with Ising model\n";
    std::cout << "3. Room assignment problem solver using simulated annealing\n";
    std::cout << "4. Exit\n";
    std::cout << "\nYour choice (1,2,3): ";
}

int main() {
    int choice;
    bool exit = false;

    while (!exit) {
        showMenu();
        std::cin >> choice;

        switch (choice) {
            case 1:
                benchmarks();
                break;
            case 2:
                isingModel();
                break;
            case 3:
                room_assignment_example();
                break;
            case 4:
                exit = true;
                std::cout << "Exit.\n";
                break;
            default:
                std::cout << "Invalid choice.\n";
                break;
        }
    }

    return 0;
}
