#include <cstddef>
#include <iostream>

#include <ctime>
#include <vector>
#include <algorithm>
#include "Data.h"

using namespace std;

#define SIZE dt->getDimension()
Data * dt = NULL;
void initializeData(char* file)
{
    cout << "Reading " << file << " ..." << endl;
    static Data data = Data(2, file);
    data.read();
    dt = &data;
    cout << "Succesfully read " << file << " !" << endl;
}

typedef struct Solution 
{
    vector<int> sequence;
    double cost;
    void print()
    {
        for(int i = 0; i < sequence.size() - 1; i++)
            cout << sequence.at(i) << " -> ";
        cout << sequence.back() << std::endl;

        cout << "Cost: " << cost << std::endl;
    }
    void print(string name) {cout<<name<<"="<<endl; print(); }
}Solution;
void updateCost(Solution &toUpdate)
{
    toUpdate.cost = 0;
    for(int i = 0; i < toUpdate.sequence.size()-1; i++) 
        toUpdate.cost += dt->getDistance(toUpdate.sequence[i], toUpdate.sequence[i+1]);
}

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
    vector<int> sorted3 = choseFromInterval(3, 2, SIZE);

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
    
    for(int i = 1; i < SIZE + 1; i++)
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
vector<InsertionInfo> calculatePossibleInsertions
    (const vector<int> & sequence, const vector<int> & inserting)
{
    std::vector<InsertionInfo> insertions = 
        std::vector<InsertionInfo>((sequence.size() - 1) * inserting.size());
        
    int count = 0;
    for(int a = 0; a < sequence.size() - 1; a++)
    {
        int predecessor = sequence[a];
        int successor = sequence[a+1];
        for(int inserted : inserting) 
        {
            insertions[count].cost = 
                - dt->getDistance(predecessor, successor)
                + dt->getDistance(predecessor, inserted)
                + dt->getDistance(inserted, successor);

            insertions[count].insertedNode = inserted;
            insertions[count].removedEdge = a;
            count++;
        }
    }
    return insertions;
}
void sortInsertions(vector<InsertionInfo> &insertions)
{
    sort(insertions.begin(), insertions.end(),
        [](InsertionInfo a, InsertionInfo b){return a.cost < b.cost;});
}
int lowerBiasedRand(int max)
{
    double alpha = (double) rand() / RAND_MAX;
    return rand() % ((int)ceil(alpha * max));
}

Solution Construct()
{
    Solution s = Solution();
    
    s.sequence = choseRandom3NodeSequence();
    
    updateCost(s);
    std::vector<int> CL = getUnusedNodes(s.sequence);
    
    while(!CL.empty())
    {
        std::vector<InsertionInfo> insertions = 
            calculatePossibleInsertions(s.sequence, CL);

        sortInsertions(insertions); 
        
        int selected = lowerBiasedRand(insertions.size());

        s.sequence.insert(
            s.sequence.begin() + insertions[selected].removedEdge + 1,
            insertions[selected].insertedNode);
        s.cost += insertions[selected].cost;

        CL = getUnusedNodes(s.sequence);
    }

    return s;
}

/*Solution ILS(int maxIter, int maxIterIls)
{
    Solution bestOfAll = Solution{{},INFINITY};
    
    for(int i = 0; i < maxIter; i++)
    {
        Solution s = Construct();

        Solution best = s;

        int iterIls = 0;

        while(iterIls <= maxIterIls)
        {
            LocalSearch(&s);
            
            if(s.cost < best.cost)
            {
                best = s;
                iterIls = 0;
            }

            s = Pertubation(best);
            
            iterIls++;
        }
        if (best.cost < bestOfAll.cost)
            bestOfAll = best;
    }
    return bestOfAll;
}*/

int main(int argc, char** argv) 
{
    srand(time(NULL));
    initializeData(argv[1]);

    cout << "Begining ILS" << endl;
    clock_t before = clock();
    Solution s = Construct();//ILS(50, SIZE/(1 + (SIZE>=150)));
    float duration = (clock()-before);
    s.print("ILS Solution");
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;


    return 0;
}