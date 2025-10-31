#include <vector>
#include <cstdlib>
#include "common.h"
#include <limits>

using namespace std;
typedef struct InsertionInfo
{
    int insertedNode;
    int removedEdge;
    double cost;
} InsertionInfo;
vector<int> choseFromInterval(int n, int first, int limit)
{
    vector<int> chosen;
    for(int i = 0; i < n; i++)
    {
        int toAdd = (rand() % ((limit-first)-i)) + first;
        bool inserted = false;
        for(int j = 0; j < chosen.size(); j++) 
        {
            if(toAdd < chosen[j]) 
            {
                chosen.insert(chosen.begin() + j,toAdd);
                inserted = true;
                break;
            }
            toAdd++;
        }
        if(!inserted) chosen.push_back(toAdd);
    }
    return chosen;
}
vector<int> choseRandom3NodeSequence()
{
    vector<int> chosen = {1};
    vector<int> sorted3 = choseFromInterval(3, 2, n);

    int i = rand()%3;
    chosen.push_back(sorted3[i]);
    sorted3.erase(sorted3.begin() + i);

    i = rand()%2;
    chosen.push_back(sorted3[i]);
    chosen.push_back(sorted3[!i]);

    chosen.push_back(1);
    
    return chosen;
}
vector<int> getUnusedNodes(const vector<int> & used)
{
    vector<int> unused;
    
    for(int i = 1; i < n + 1; i++)
    {
        bool isUsed = false;
        for(int node : used)
        {
            if(i == node) 
            {
                isUsed = true;
                break;
            }   
        }
        if(!isUsed) unused.push_back(i);
    }

    return unused;
}
InsertionInfo calculateBestInsertion
    (const vector<int> & sequence, const vector<int> & inserting)
{
    InsertionInfo best;
    best.cost = numeric_limits<double>::infinity();

    for(int a = 0; a < sequence.size() - 1; a++)
    {
        int predecessor = sequence[a];
        int successor = sequence[a+1];
        for(int inserted : inserting) 
        {
			double cost = 
                - getMatrixCost(predecessor-1, successor-1)
                + getMatrixCost(predecessor-1, inserted-1)
                + getMatrixCost(inserted-1, successor-1);

			if(cost < best.cost)
			{
				best.cost = cost;
				best.insertedNode = inserted;
				best.removedEdge = a;
			}
        }
    }
    return best;
}

Solution Construct()
{
    Solution s = Solution();
    
    s.sequence = choseRandom3NodeSequence();
    
    updateCost(s);
    std::vector<int> CL = getUnusedNodes(s.sequence);
    
    while(!CL.empty())
    {
        InsertionInfo insertion = 
            calculateBestInsertion(s.sequence, CL);
        
        s.sequence.insert(
            s.sequence.begin() + insertion.removedEdge + 1,
            insertion.insertedNode);
        s.cost += insertion.cost;

        CL = getUnusedNodes(s.sequence);
    }

    return s;
}
