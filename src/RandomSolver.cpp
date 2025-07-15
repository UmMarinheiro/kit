#pragma once

#include "Data.h"
#include "Solver.h"

class RandomSolver : public Solver
{
public:
    RandomSolver(Data *_data) : Solver(_data)
    {
        solution.sequence.push_back(1);
        for(int i = 0; i < dimension-1; i++)
        {
            vector<int> unused = getUnusedNodes(solution.sequence);
            solution.sequence.push_back(unused.at(rand()%unused.size()));
        }
        solution.sequence.push_back(1);

        updateCost();
    }
};