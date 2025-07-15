#pragma once

#include "Solver.h"
#include <vector>
#include <algorithm>

typedef struct InsertionInfo
{
    int insertedNode;
    int removedEdge;
    double cost;
} InsertionInfo;

class ILSSolver : public Solver
{
public:
    ILSSolver(Data *_data, int maxIter, int maxIterIls);
private:

    vector<int> choseFromInterval(int n, int first, int size);
    vector<int> choseRandom3NodeSequence();
    std::vector<InsertionInfo> calculatePossibleInsertions
        (const vector<int> & sequence, const vector<int> & inserting);
    void sortInsertions(vector<InsertionInfo> &insertions);
    int lowerBiasedRand(int max);
    
    Solution Construct();


    bool bestImprovementSwap(Solution *s);
    bool bestImprovement2Opt(Solution *s);
    bool bestImprovementOrOpt(Solution *s, int nVertex);
    void LocalSearch(Solution *s);
        
    Solution Pertubation(const Solution &s);
};
