#ifndef ROOM_ASSIGNMENT_H
#define ROOM_ASSIGNMENT_H

#include <vector>
#include <random>
#include <utility>

class RoomAssignment {
public:
    // Constructor
    RoomAssignment(int num_students, const std::vector<std::vector<double>>& dislike_matrix);

    // Main solver method
    void solve();

    // Utility methods
    double calculateTotalDislikes();
    void printAssignments();

private:
    // Member variables
    int n;  // number of students
    std::vector<std::vector<double>> dislikes;  // dislike matrix
    std::vector<int> assignments;  // current room assignments
    std::random_device rd;
    std::mt19937 gen;

    // Helper methods
    std::pair<int, int> selectRandomStudents();
};

void room_assignment_example();

#endif // ROOM_ASSIGNMENT_H
