#include <algorithm>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include "Data.h"

using namespace std;

typedef struct Solution
{
    std::vector<int> sequence;
    double cost;

    Solution() : sequence({}), cost(0.0){}
    Solution(vector<int> s,double c) : sequence(s), cost(c){}

    const void print()
    {
        for(int i = 0; i < sequence.size() - 1; i++)
        {
            if(sequence.at(i)<100)
            {
                std::cout << " ";
                if(sequence.at(i)<10) std::cout << " ";
            }
            
            std::cout << sequence.at(i) << " -> ";
        }
        std::cout << sequence.back() << std::endl;

        cout << "Cost: " << cost << std::endl;
    }
    const void print(char* name)
    {
        cout << name << " = " << std::endl;
        print();
    }
}Solution;

typedef struct InsertionInfo
{
    int insertedNode;
    int removedEdge;
    double cost;
} InsertionInfo;

class Solver
{
    public:
        Solver(Data *_data) : data(_data){}
        
        Solution solution;

    protected: 
        void updateCost()
        {
            updateCost(solution);
        }
        void updateCost(Solution &s)
        {
            s.cost = 0;
            auto matrixCost = data->getMatrixCost();
            for(int i = 0; i < s.sequence.size() - 1; i++)
            {
                s.cost += matrixCost[s.sequence[i]-1][s.sequence[i+1]-1];
            }
        }
        Data *data = NULL;
};

class CrescentSolver : public Solver
{
    public:
        CrescentSolver(Data *_data) : Solver(_data)
        {
            for(int i = 1; i < data->getDimension()+1; i++)
                solution.sequence.push_back(i);
            solution.sequence.push_back(1);

            updateCost();
        }
};
class RandomSolver : public Solver
{
    public:
        RandomSolver(Data *_data) : Solver(_data)
        {
            solution.sequence.push_back(1);
            for(int i = 1; i < data->getDimension(); i++)
            {
                vector<int> unused = getUnusedNodes(solution.sequence);
                solution.sequence.push_back(unused.at(rand()%unused.size()));
            }
            solution.sequence.push_back(1);

            updateCost();
        }
    private:
        vector<int> getUnusedNodes(const vector<int> & used)
        {
            vector<int> unused;
            
            for(int i = 1; i <= data->getDimension(); i++)
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
};

class ILSSolver : public Solver
{
    public:
        ILSSolver(Data *_data, int maxIter, int maxIterIls) : Solver(_data)
        {
            Solution bestOfAll = Solution({},INFINITY);
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
            solution = bestOfAll;
        }
        ILSSolver(Data *_data) : Solver(_data)
        {
            cout<<"Running ILS..."<<endl;
            solution = Construct();
            solution.print((char*)"ILS Pre-Search");
            LocalSearch(&solution);
        }
    private:
        vector<int> choseRandom3NodeSolution()
        {
            vector<int> chosen = 
            {
                1,
                rand()%(data->getDimension() - 1) + 2,
                rand()%(data->getDimension() - 2) + 2,
                rand()%(data->getDimension() - 3) + 2,
                1
            };

            if(chosen[2]>=chosen[1]) chosen[1]++;

            if(chosen[3]>=chosen[1]) chosen[2]++;
            if(chosen[3]>=chosen[2]) chosen[2]++;
            
            return chosen;
        }
        vector<int> getUnusedNodes(const vector<int> & used)
        {
            vector<int> unused;
            
            for(int i = 1; i <= data->getDimension(); i++)
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
        std::vector<InsertionInfo> calculatePossibleInsertions
            (const vector<int> & sequence, const vector<int> & inserting)
        {
            auto adjMatrix = data->getMatrixCost();
             
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
                        - adjMatrix[predecessor - 1][successor - 1]
                        + adjMatrix[successor - 1][inserted - 1]
                        + adjMatrix[predecessor - 1][inserted - 1];

                    insertions[count].insertedNode = inserted;
                    insertions[count].removedEdge = a;
                    count++;
                }
            }
            return insertions;
        }
        void sortInsertions(vector<InsertionInfo> &insertions)
        {
            sort(insertions.begin(),insertions.end(),
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

            s.sequence = choseRandom3NodeSolution();
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
                CL = getUnusedNodes(s.sequence);
            }
            
            updateCost(s);
            return s;
        }

        bool bestImprovementSwap(Solution *s)
        {
            auto adjMatrix = data->getMatrixCost();

            double bestDelta = 0;
            int best_i, best_j;

            for(int i = 1; i < s->sequence.size() -1; i++)
            {
                int vi = s->sequence[i];
                int vi_predecessor = s->sequence[i-1];
                int vi_succesor = s->sequence[i+1];

                for(int j = i+1; j < s->sequence.size() - 1; j++)
                {
                    int vj = s->sequence[j];
                    int vj_predecessor = s->sequence[j-1];
                    int vj_succesor = s->sequence[j+1];

                    double delta = 
                        -adjMatrix[vi_predecessor-1][vi-1]
                        -adjMatrix[vi-1][vi_succesor-1] 
                        +adjMatrix[vi_predecessor-1][vj-1]
                        +adjMatrix[vj-1][vi_succesor-1]
                        
                        -adjMatrix[vj_predecessor-1][vj-1] 
                        -adjMatrix[vj-1][vj_succesor-1]
                        +adjMatrix[vj_predecessor-1][vi-1] 
                        +adjMatrix[vi-1][vj_succesor-1];

                    if(vi == vj_predecessor) 
                        delta -= 
                            -adjMatrix[vi-1][vi_succesor-1]
                            -adjMatrix[vj_predecessor-1][vj-1];

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
            auto adjMatrix = data->getMatrixCost();

            double bestDelta = 0;
            int best_i, best_j;

            for(int i = 1; i < s->sequence.size(); i++)
            {
                int ei_start = s->sequence[i-1];
                int ei_end = s->sequence[i];

                for(int j = i+2; j < s->sequence.size(); j++)
                {
                    int ej_start = s->sequence[j-1];
                    int ej_end = s->sequence[j];

                    double delta = 
                        -adjMatrix[ei_start-1][ei_end-1]
                        -adjMatrix[ej_start-1][ej_end-1]
                        
                        +adjMatrix[ei_start-1][ej_start-1]
                        +adjMatrix[ei_end-1][ej_end-1];

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
                for(int i = best_i, j = best_j-1; i < j; i++, j--)
                {
                    std::swap(s->sequence[i],s->sequence[j]);
                }
                s->cost = s->cost + bestDelta;

                return true;
            }
            else return false;
        }
        bool bestImprovementOrOpt(Solution *s, int nVertex)
        {
            auto adjMatrix = data->getMatrixCost();

            double bestDelta = 0;
            int best_i, best_j;

            for(int i = 1; i < s->sequence.size() - 1 - (nVertex-1); i++)
            {
                int block_predecessor = s->sequence[i-1];
                int block_start = s->sequence[i];
                int block_end = s->sequence[i+(nVertex-1)];
                int block_succesor = s->sequence[i+(nVertex-1) +1];

                for(int j = 1; j < s->sequence.size(); j++)
                {
                    int edge_start = s->sequence[j-1];
                    int edge_end = s->sequence[j];

                    if(i-1 <= j-1 && j <= i+(nVertex-1) +1) continue;
                    // block_predecessor_index <= edge_start_index && 
                    // edge_end_index <= block_succsor_index 
                    // ou seja, aresta não faz contato com o bloco

                    double delta = 
                        -adjMatrix[block_predecessor-1][block_start-1] 
                        -adjMatrix[block_end-1][block_succesor-1]
                        
                        +adjMatrix[block_predecessor-1][block_succesor-1]

                        -adjMatrix[edge_start-1][edge_end-1]
                        +adjMatrix[edge_start-1][block_start-1]
                        +adjMatrix[block_end-1][edge_end-1];


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
                for(int i = 0; i < nVertex;i++)
                {
                    s->sequence.insert(s->sequence.begin() + best_j  + (i)*(best_j<best_i),
                        s->sequence[best_i + (i)*(best_j<best_i)]);
                    s->sequence.erase(s->sequence.begin() + best_i + 
                        (1+i)*(best_j<best_i));
                }
                s->cost = s->cost + bestDelta;

                return  true;
            }
            else return false;
        }
        void LocalSearch(Solution *s)
        {
            std::vector<int> NL = {1, 2, 3, 4,5};
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
                
                
                if(improved) NL = {1,2,3,4,5};
                else NL.erase(NL.begin() + n);
            }
        } //TODO
        Solution Pertubation(const Solution &s){return s;} //TODO
};

int main(int argc, char** argv) 
{
    srand(time(0));

    cout << "Reading " << argv[1] << " ..." << endl;

    Data data = Data(argc, argv[1]);
    data.read();

    cout << "Succesfully read " << argv[1] << " !" << endl;

    Solver *solver = new ILSSolver(&data);
    solver->solution.print((char*)"ILSSolution");
    delete(solver);

    solver = new RandomSolver(&data);
    solver->solution.print((char*)"RandomSolution");
    delete(solver);

    return 0;
}