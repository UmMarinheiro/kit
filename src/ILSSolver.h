#pragma once

#include "Solution.h"
#include "Solver.h"
#include <vector>
#include <algorithm>

typedef struct Subsequence
{
    double T, C;
    int W;
    int first, last;
} Subsequence;
typedef struct BuildedSolution
{
    Solution solution;
    vector<vector<Subsequence>> subseq_matrix;
    Subsequence getSubsequence(int startingNode, int endingNode) const;
} BuildedSolution;

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
    inline Subsequence Concatenate(Subsequence s1, Subsequence s2);
    void UpdateAllSubseq(BuildedSolution &b);

    vector<int> choseFromInterval(int n, int first, int size);
    vector<int> choseRandom3NodeSequence();
    std::vector<InsertionInfo> calculatePossibleInsertions
        (const vector<int> & sequence, const vector<int> & inserting);
    void sortInsertions(vector<InsertionInfo> &insertions);
    int lowerBiasedRand(int max);
    
    BuildedSolution Construct();


    bool bestImprovementSwap(BuildedSolution *s);
    bool bestImprovement2Opt(BuildedSolution *s);
    bool bestImprovementOrOpt(BuildedSolution *s, int nVertex);
    void LocalSearch(BuildedSolution *s);
        
    BuildedSolution Pertubation(const BuildedSolution &s);
};
