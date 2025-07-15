#pragma once

#include "Data.h"
#include "Solution.h"
#include <vector>

#define FIRST_NODE 1

class Solver
{
    public:
        Solver(Data *_data);
        Solution solution;

    protected: 
        Data *data;
        int dimension;
        
        double getCost(int node1, int node2);
        void updateCost();
        void updateCost(Solution &s);
        vector<int> getUnusedNodes(const vector<int> &used);
};