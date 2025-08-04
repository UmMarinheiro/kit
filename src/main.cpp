#include <iostream>

#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
#include "Data.h"

using namespace std;

#define SIZE dt->getDimension()
Data * dt = NULL;
void initializeData(char* file)
{
    cout << "Reading " << file << " ..." << endl;
    dt = new Data(2, file);
    dt->read();
    cout << "Succesfully read " << file << " !" << endl;
}
void freeData() {delete dt;}

typedef struct Solution 
{
    vector<int> sequence;
    double cost;
    void print() const
    {
        for(int i = 0; i < sequence.size() - 1; i++)
            cout << sequence.at(i) << " -> ";
        cout << sequence.back() << std::endl;

        cout << "Cost: " << cost << std::endl;
    }
    void print(string name) const {cout<<name<<"="<<endl; print(); }
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
    int limitOfBest = (int)ceil(alpha * max);
    return limitOfBest>0?rand()%limitOfBest : 0;
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

bool bestImprovementSwap(Solution *s)
{
    double bestDelta = 0;
    int best_i, best_j;

    for(int i = 1; i < s->sequence.size() - 1; i++)
    {
        int vi = s->sequence[i];
        int vi_predecessor = s->sequence[i-1];
        int vi_succesor = s->sequence[i+1];

        for(int j = i+2; j < s->sequence.size() - 1; j++)
        {
            int vj = s->sequence[j];
            int vj_predecessor = s->sequence[j-1];
            int vj_succesor = s->sequence[j+1];

            double delta = 
                -dt->getDistance(vi_predecessor, vi)
                -dt->getDistance(vi, vi_succesor) 
                +dt->getDistance(vi_predecessor, vj)
                +dt->getDistance(vj, vi_succesor)
                
                -dt->getDistance(vj_predecessor, vj)
                -dt->getDistance(vj, vj_succesor)
                +dt->getDistance(vj_predecessor, vi)
                +dt->getDistance(vi, vj_succesor);

            if(delta < bestDelta)
            {
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    if(bestDelta < 0) 
    {
        std::swap(s->sequence[best_i], s->sequence[best_j]);
        s->cost = s->cost + bestDelta;

        return  true;
    }
    else return false;
}
bool bestImprovement2Opt(Solution *s)
{
    double bestDelta = 0;
    int best_i, best_j;

    for(int i = 0; i < s->sequence.size()-1; i++)
    {
        int ei_start = s->sequence[i];
        int ei_end = s->sequence[i+1];

        for(int j = i+2; j < s->sequence.size()-1; j++)
        {
            int ej_start = s->sequence[j];
            int ej_end = s->sequence[j+1];

            double delta = 
                -dt->getDistance(ei_start, ei_end)
                -dt->getDistance(ej_start, ej_end)
                
                +dt->getDistance(ei_start, ej_start)
                +dt->getDistance(ei_end, ej_end);

            if(delta < bestDelta)
            {
                bestDelta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(bestDelta < 0) 
    {
        reverse(s->sequence.begin() + best_i + 1, s->sequence.begin() + best_j + 1);
        s->cost = s->cost + bestDelta;

        return true;
    }
    else return false;
}
bool bestImprovementOrOpt(Solution *s, int nVertex)
{
    double bestDelta = 0;
    int best_block, best_edge;

    for(int block_start_index = 1; block_start_index < s->sequence.size()-1 - (nVertex-1); block_start_index++)
    {
        int block_predecessor = s->sequence[block_start_index-1];
        int block_start = s->sequence[block_start_index];
        int block_end = s->sequence[block_start_index + nVertex-1];
        int block_succesor = s->sequence[block_start_index + nVertex];
        
        for(int j = 0; j < s->sequence.size()-1 - (nVertex+2); j++)
        {
            int edge_start_index = j + (j >= block_start_index-1)*(nVertex+2);

            int edge_start = s->sequence[edge_start_index];
            int edge_end = s->sequence[edge_start_index+1];

            double delta = 
                -dt->getDistance(block_predecessor, block_start) 
                -dt->getDistance(block_end, block_succesor)
                
                +dt->getDistance(block_predecessor, block_succesor)
                -dt->getDistance(edge_start, edge_end)
                
                +dt->getDistance(edge_start, block_start)
                +dt->getDistance(block_end, edge_end);


            if(delta < bestDelta)
            {
                bestDelta = delta;
                best_block = block_start_index;
                best_edge = edge_start_index;
            }
        }
    }

    if(bestDelta < 0) 
    {
        if(best_block<best_edge)
            rotate(
                s->sequence.begin()+best_block, 
                s->sequence.begin()+best_block + nVertex, 
                s->sequence.begin()+best_edge + 1);
        else 
            rotate(
                s->sequence.begin()+best_edge + 1, 
                s->sequence.begin()+best_block, 
                s->sequence.begin()+best_block + nVertex);
        s->cost = s->cost + bestDelta;

        return  true;
    }
    else return false;
}

void LocalSearch(Solution *s)
{
    std::vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while(!NL.empty())
    {
        int n = rand() % NL.size();

        switch (NL[n]) 
        {
        case 1:
            improved = bestImprovementSwap(s);
            break;
        case 2:
            improved = bestImprovement2Opt(s);
            break;
        case 3:
            improved = bestImprovementOrOpt(s,1);
            break;
        case 4:
            improved = bestImprovementOrOpt(s,2);
            break;
        case 5:
            improved = bestImprovementOrOpt(s,3);
            break;
        }

        if(improved) 
        {
            NL = {1,2,3,4,5};
        }
        else NL.erase(NL.begin() + n);
    }
} 
Solution Pertubation(const Solution &s)
{
    vector<int> delimiters = choseFromInterval(4, 1, SIZE);
    
    Solution sr = {{}, 0};
    sr.sequence.insert(sr.sequence.end(), s.sequence.begin(), s.sequence.begin() + delimiters[0]);
    sr.sequence.insert(sr.sequence.end(), s.sequence.begin() + delimiters[2], s.sequence.begin() + delimiters[3]);
    sr.sequence.insert(sr.sequence.end(), s.sequence.begin() + delimiters[1], s.sequence.begin() + delimiters[2]);
    sr.sequence.insert(sr.sequence.end(), s.sequence.begin() + delimiters[0], s.sequence.begin() + delimiters[1]);
    sr.sequence.insert(sr.sequence.end(), s.sequence.begin() + delimiters[3], s.sequence.end());
    
    double delta = 
        -dt->getDistance(s.sequence[delimiters[0]-1], s.sequence[delimiters[0]])
        -dt->getDistance(s.sequence[delimiters[1]-1], s.sequence[delimiters[1]])
        +dt->getDistance(s.sequence[delimiters[2]-1], s.sequence[delimiters[0]])
        +dt->getDistance(s.sequence[delimiters[1]-1], s.sequence[delimiters[3]])
        
        -dt->getDistance(s.sequence[delimiters[2]-1], s.sequence[delimiters[2]])
        -dt->getDistance(s.sequence[delimiters[3]-1], s.sequence[delimiters[3]])
        +dt->getDistance(s.sequence[delimiters[0]-1], s.sequence[delimiters[2]])
        +dt->getDistance(s.sequence[delimiters[3]-1], s.sequence[delimiters[1]]);

    sr.cost = s.cost + delta;

    return sr;
}

Solution ILS(int maxIter, int maxIterIls)
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
}

int main(int argc, char** argv) 
{
    srand(time(NULL));
    
	initializeData(argv[1]);
	
	cout << "Running ILS..." << endl;
	clock_t before = clock();
	Solution s = ILS(50, SIZE/(1 + (SIZE>=150)));            
	
	float duration = (clock()-before);
	s.print("ILS Solution");
	cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;
	
    freeData();

    return 0;
}