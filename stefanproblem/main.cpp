#include "stefanproblem/stefanProblem.h"
#include <iostream>

int main(int argc, char *argv[])
{
    std::cout << " Start calculations...\n\n";
    StefanProblemSolver solver;
    solver.calculate(21);

    std::cout << " Solve for 21 pts:\n";
    for (int i = 0; i < solver.PTB.size(); ++i) {
        std::cout << "ts - " << solver.ts[i] << " ms: " << solver.PTB[i] << "(PTB)\n";
    }    

    solver.clean();
    solver.calculate(51);

    std::cout << " Solve for 51 pts:\n";
    for (int i = 0; i < solver.PTB.size(); ++i) {
        std::cout << "ts - " << solver.ts[i] << " ms: " << solver.PTB[i] << "(PTB)\n";
    }

    /*solver.clean();
    solver.calculate(101);

    std::cout << " Solve for 101 pts:\n";
    for (int i = 0; i < solver.PTB.size(); ++i) {
        std::cout << "ts - " << solver.ts[i] << " ms: " << solver.PTB[i] << "(PTB)\n";
    }

    solver.clean();
    solver.calculate(201);

    std::cout << " Solve for 201 pts:\n";
    for (int i = 0; i < solver.PTB.size(); ++i) {
        std::cout << "ts - " << solver.ts[i] << " ms: " << solver.PTB[i] << "(PTB)\n";
    }*/

    return 0;
}
