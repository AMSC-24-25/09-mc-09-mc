#include "benchmarks.h"
#include "ising_model.h"

#include <iostream>

void showMenu() {
    std::cout << "\nWhat would you like to run?\n";
    std::cout << "1. Benchmarks for each integration and domain\n";
    std::cout << "2. Specific heat per particle with Ising model\n";
    std::cout << "3. Exit\n";
    std::cout << "Your choice (1,2,3): ";
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