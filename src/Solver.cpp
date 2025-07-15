#include "Data.h"
#include "Solution.h"

#include "Solver.h"
    

Solver::Solver(Data *_data) : data(_data), dimension(_data->getDimension()){};
double Solver::getCost(int node1, int node2)
{
    return data->getMatrixCost()[node1 - FIRST_NODE][node2 - FIRST_NODE];
}
void Solver::updateCost()
{
    updateCost(solution);
}
void Solver::updateCost(Solution &s)
{
    s.cost = 0;
    for(int i = 0; i < s.sequence.size() - 1; i++)
    {
        s.cost += getCost(s.sequence[i], s.sequence[i+1]);
    }
}
vector<int> Solver::getUnusedNodes(const vector<int> & used)
{
    vector<int> unused;
    
    for(int i = FIRST_NODE; i < dimension + FIRST_NODE; i++)
    {
        bool isUsed = false;
        for(int j = 0; j < used.size(); j++)
        {
            if(i == used[j]) 
            {
                isUsed = true;
                break;
            }   
        }
        if(!isUsed) unused.push_back(i);
    }

    return unused;
}