#include <iostream>

#include <ctime>

#include "Solution.h"
#include "Solver.h"
#include "RandomSolver.cpp"
#include "ILSSolver.h"
#include "Data.h"

using namespace std;

int main(int argc, char** argv) 
{
    srand(time(NULL));

    cout << "Reading " << argv[1] << " ..." << endl;

    Data data = Data(argc, argv[1]);
    data.read();

    cout << "Succesfully read " << argv[1] << " !" << endl;
    
    cout << "Begining ILS" << endl;
    clock_t before = clock();
    Solver *solver = new ILSSolver(&data,
        50,
        data.getDimension()/(1 + (data.getDimension()>=150)));
    float duration = (clock()-before);
    solver->solution.print((char*)"ILSSolution");
    cout << "Concluded ILS - Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;
    delete(solver);

    /*solver = new RandomSolver(&data);
    solver->solution.print((char*)"RandomSolution");
    delete(solver);*/

    return 0;
}