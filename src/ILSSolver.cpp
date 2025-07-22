#include "ILSSolver.h"

#include "Solution.h"
#include "Solver.h"
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>

#include <cfloat>

Subsequence BuildedSolution::getSubsequence(int startingNode, int endingNode) const {return subseq_matrix[startingNode-FIRST_NODE][endingNode-FIRST_NODE];}

inline Subsequence ILSSolver::Concatenate(Subsequence s1, Subsequence s2)
{
    Subsequence s;
    double temp = getCost(s1.last, s2.first);
    s.W = s1.W + s2.W;
    s.T = s1.T + temp + s2.T;
    s.C = s1.C + s2.W*(s1.T + temp) + s2.C;
    s.first = s1.first;
    s.last = s2.last;

    return  s;
}
void ILSSolver::UpdateAllSubseq(BuildedSolution &b)
{
    int n = b.solution.sequence.size();
    b.subseq_matrix = vector<vector<Subsequence>>(n, vector<Subsequence>(n));

    for(int i = 0; i < n; i++)
    {
        b.subseq_matrix[i][i].W = (i>0);
        b.subseq_matrix[i][i].C = 0;
        b.subseq_matrix[i][i].T = 0;

        b.subseq_matrix[i][i].first = b.solution.sequence[i];
        b.subseq_matrix[i][i].last = b.solution.sequence[i];
    }

    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
            b.subseq_matrix[i][j] = Concatenate(b.subseq_matrix[i][j-1],b.subseq_matrix[j][j]);

    for(int i = n-1; i >= 0; i--)
        for(int j = i-1; j >=0; j--)
            b.subseq_matrix[i][j] = Concatenate(b.subseq_matrix[i][j+1],b.subseq_matrix[j][j]);

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            std::cout << b.subseq_matrix[i][j].C << " ";
        std::cout<<"\n";
    }
    updateCost(b.solution);
    double cost = 0;
    std::cout<<"\n\n\n";
    for(int i = 0; i < dimension-1; i++)
    {
        std::cout<<cost<<" "<<b.getSubsequence(1, i+1).C<<"\n";
        cost += cost + getCost(b.solution.sequence[i], b.solution.sequence[i+1]);
    }
    getchar();
}

ILSSolver::ILSSolver(Data *_data, int maxIter, int maxIterIls) : Solver(_data)
{
    Solution bestOfAll = Solution({},INFINITY);
    
    for(int i = 0; i < maxIter; i++)
    {
        std::cout<<"Starting construction\n";
        BuildedSolution s = Construct();
        std::cout<<"Concluded construction\n";

        BuildedSolution best = s;

        int iterIls = 0;

        while(iterIls <= maxIterIls)
        {

            std::cout<<"Starting LS\n";
            LocalSearch(&s);
            std::cout<<"Concluded LS\n";
            
            if(s.solution.cost < best.solution.cost)
            {
                best = s;
                iterIls = 0;
            }

            std::cout<<"Starting Pertubation\n";
            s = Pertubation(best);
            std::cout<<"Concluded Pertubation\n";
            
            iterIls++;
        }
        if (best.solution.cost < bestOfAll.cost)
            bestOfAll = best.solution;
    }
    solution = bestOfAll;
    solution.print((char*)"Antes");
    updateCost(solution);
    solution.print((char*)"Depois");
}
vector<int> ILSSolver::choseFromInterval(int n, int first, int size)
{
    vector<int> chosen;
    for(int i = first; i < n + first; i++)
    {
        int toAdd = (rand() % (size-(i-first))) + first;
        char inserted = 0;
        for(int j = 0; j < chosen.size(); j++) 
        {
            if(toAdd < chosen[j]) 
            {
                chosen.insert(chosen.begin() + j,toAdd);
                inserted = 1;
                break;
            }
            toAdd++;
        }
        if(!inserted) chosen.push_back(toAdd);
    }
    return chosen;
}
vector<int> ILSSolver::choseRandom3NodeSequence()
{
    vector<int> chosen = {FIRST_NODE};
    vector<int> sorted3 = choseFromInterval(3,FIRST_NODE + 1,dimension - 1);

    int i = rand()%3;
    chosen.push_back(sorted3[i]);
    sorted3.erase(sorted3.begin() + i);

    i = rand()%2;
    chosen.push_back(sorted3[i]);
    sorted3.erase(sorted3.begin() + i);

    chosen.push_back(sorted3[0]);
    chosen.push_back(1);
    
    return chosen;
}
std::vector<InsertionInfo> ILSSolver::calculatePossibleInsertions
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
                - getCost(predecessor, successor)
                + getCost(predecessor, inserted)
                + getCost(inserted, successor);

            insertions[count].insertedNode = inserted;
            insertions[count].removedEdge = a;
            count++;
        }
    }
    return insertions;
}
void ILSSolver::sortInsertions(vector<InsertionInfo> &insertions)
{
    sort(insertions.begin(),insertions.end(),
        [](InsertionInfo a, InsertionInfo b){return a.cost < b.cost;});
}
int ILSSolver::lowerBiasedRand(int max)
{
    double alpha = (double) rand() / RAND_MAX;
    return rand() % ((int)ceil(alpha * max));
}

BuildedSolution ILSSolver::Construct()
{
    BuildedSolution b = BuildedSolution();
    Solution & s = b.solution;
    
    s.sequence = choseRandom3NodeSequence();
    
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

    UpdateAllSubseq(b);
    s.cost = b.getSubsequence(1, dimension).C;
    return b;
}

bool ILSSolver::bestImprovementSwap(BuildedSolution *b)
{
    std::cout<<"Starting SWAP\n";
    Solution *s = &(b->solution);
    double bestCost = DBL_MAX;
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

            Subsequence sb = b->getSubsequence(1, vi_predecessor);
            sb = Concatenate(sb, b->getSubsequence(vj, vj));
            sb = Concatenate(sb, b->getSubsequence(vi_succesor, vj_predecessor));
            sb = Concatenate(sb, b->getSubsequence(vi, vi));
            sb = Concatenate(sb, b->getSubsequence(vj_succesor, dimension));
            
            if(sb.C < b->solution.cost)
            {
                bestCost = sb.C;
                best_i = i;
                best_j = j;
            }
        }
    }
    std::cout<<"Ended SWAP\n";
    
    if(bestCost < b->solution.cost) 
    {
        std::swap(s->sequence[best_i], s->sequence[best_j]);
        s->cost = bestCost;

        cout<<"Trocando "<<s->sequence[best_i] << " e " << s->sequence[best_j]<<"\n";
        UpdateAllSubseq(*b);

        return  true;
    }
    else return false;
}
bool ILSSolver::bestImprovement2Opt(BuildedSolution *b)
{
    std::cout<<"Starting 2Opt\n";
    Solution *s = &(b->solution);
    double bestCost = DBL_MAX;
    int best_i, best_j;

    for(int i = 0; i < s->sequence.size()-1; i++)
    {
        int ei_start = s->sequence[i];
        int ei_end = s->sequence[i+1];

        for(int j = i+2; j < s->sequence.size()-1; j++)
        {
            int ej_start = s->sequence[j];
            int ej_end = s->sequence[j+1];

            Subsequence sb = b->getSubsequence(1, ei_start);
            sb = Concatenate(sb, b->getSubsequence(ej_start, ei_end));
            sb = Concatenate(sb, b->getSubsequence(ei_end, dimension));
            
            if(sb.C < bestCost)
            {
                bestCost = sb.C;
                best_i = i;
                best_j = j;
            }
        }
    }

    std::cout<<"Ended 2OPt\n";

    if(bestCost < 0) 
    {
        for(int i = best_i+1, j = best_j; i < j; i++, j--)
        {
            std::swap(s->sequence[i],s->sequence[j]);
        }
        s->cost = bestCost;
        cout<<"Invertendo de "<<s->sequence[best_i+1] << " a " << s->sequence[best_j] <<"\n";
        UpdateAllSubseq(*b);
        
        return true;
    }
    else return false;
}
bool ILSSolver::bestImprovementOrOpt(BuildedSolution *b, int nVertex)
{
    Solution *s = &(b->solution);
    double bestCost = 0;
    int best_block, best_edge;

    for(int block_start_index = 1; block_start_index < s->sequence.size()-1 - (nVertex-1); block_start_index++)
    {
        int block_predecessor = s->sequence[block_start_index-1];
        int block_start = s->sequence[block_start_index];
        int block_end = s->sequence[block_start_index+(nVertex-1)];
        int block_succesor = s->sequence[block_start_index+(nVertex-1) +1];
        
        for(int j = 0; j < s->sequence.size()-1 - (nVertex+2); j++)
        {
            int edge_start_index = j + (j >= block_start_index-1)*(nVertex+2);

            int edge_start = s->sequence[edge_start_index];
            int edge_end = s->sequence[edge_start_index+1];

            Subsequence sb;
            if(j >= block_start_index-1)
            {
                sb = b->getSubsequence(1, block_predecessor);
                sb = Concatenate(sb, b->getSubsequence(block_succesor, edge_start));
                sb = Concatenate(sb, b->getSubsequence(block_start, block_end));
                sb = Concatenate(sb, b->getSubsequence(edge_end, dimension));
            }
            else 
            {
                sb = b->getSubsequence(1, edge_start);
                sb = Concatenate(sb, b->getSubsequence(block_start, block_end));
                sb = Concatenate(sb, b->getSubsequence(edge_end, block_start));
                sb = Concatenate(sb, b->getSubsequence(block_end, dimension));
            }
            if(sb.C < bestCost)
            {
                bestCost = sb.C;
                best_block = block_start_index;
                best_edge = edge_start_index;
            }
        }
    }

    if(bestCost < 0) 
    {
        for(int i = 0; i < nVertex;i++)
        {       
            s->sequence.insert(s->sequence.begin() + (best_edge+1) + (i)*(best_edge<best_block),
                s->sequence[best_block + (i)*(best_edge<best_block)]);
            s->sequence.erase(s->sequence.begin() + best_block + 
                (1+i)*(best_edge<best_block));
        }
        s->cost = bestCost;
        cout<<"Movendo de "<<s->sequence[best_block] << " a " << s->sequence[best_block + nVertex - 1] << " para depois de " << s->sequence[best_edge] <<"\n";
        UpdateAllSubseq(*b);

        return  true;
    }
    else return false;
}
void ILSSolver::LocalSearch(BuildedSolution *b)
{
    std::vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while(!NL.empty())
    {
        int n = rand() % NL.size();

        switch (NL[n]) 
        {
        case 1:
            improved = bestImprovementSwap(b);
            break;
        case 2:
            improved = bestImprovement2Opt(b);
            break;
        case 3:
            improved = bestImprovementOrOpt(b,1);
            break;
        case 4:
            improved = bestImprovementOrOpt(b,2);
            break;
        case 5:
            improved = bestImprovementOrOpt(b,3);
            break;
        }

        if(improved) 
        {
            NL = {1,2,3,4,5};
        }
        else NL.erase(NL.begin() + n);
    }
} 
BuildedSolution ILSSolver::Pertubation(const BuildedSolution &b)
{
    const Solution &s = b.solution;
    vector<int> delimiters = choseFromInterval(4, 1, s.sequence.size()-1);

    for(int i = 0; i < delimiters.size(); i++) cout<<delimiters[i] <<" ";
    cout<<"\n";

    BuildedSolution br;
    Solution &sr = br.solution;

    for(int i = 0; i < delimiters[0]; i++) sr.sequence.push_back(s.sequence.at(i));
    for(int i = delimiters[2]; i < delimiters[3]; i++) sr.sequence.push_back(s.sequence.at(i));
    for(int i = delimiters[1]; i < delimiters[2]; i++) sr.sequence.push_back(s.sequence.at(i));
    for(int i = delimiters[0]; i < delimiters[1]; i++) sr.sequence.push_back(s.sequence.at(i));
    for(int i = delimiters[3]; i < s.sequence.size(); i++) sr.sequence.push_back(s.sequence.at(i));
    
    Subsequence sb = b.getSubsequence(1,(delimiters[0]-1 > 1? delimiters[0]-1 : 1));
    sb = Concatenate(sb, b.getSubsequence(delimiters[2],delimiters[3]-1));
    sb = Concatenate(sb, b.getSubsequence(delimiters[1],delimiters[2]-1));
    sb = Concatenate(sb, b.getSubsequence(delimiters[0],delimiters[1]-1));
    sb = Concatenate(sb, b.getSubsequence(delimiters[3],dimension));

    sr.cost = sb.C;
    UpdateAllSubseq(br);

    return br;
}
