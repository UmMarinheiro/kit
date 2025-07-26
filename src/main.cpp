#include <cstdio>
#include <iostream>

#include <iomanip>
#include <ctime>
#include <ostream>
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
    double duration = 0;
    toUpdate.cost = 0;
    
    for(int i = 0; i < toUpdate.sequence.size()-1; i++)
    {
        duration += dt->getDistance(toUpdate.sequence[i], toUpdate.sequence[i+1]);
        toUpdate.cost += duration;
    }
}
void verifyCost(Solution &toVerify)
{
    Solution s(toVerify);
    updateCost(s);
    cout<<toVerify.cost-s.cost<<endl;
}

typedef struct Subsequence
{
    double t,c;
    int w;
    int first, last;
    inline static Subsequence Concatenate(Subsequence &s1, Subsequence &s2)
    {
        Subsequence s;
        double temp = dt->getDistance(s1.last, s2.first);
        s.w = s1.w + s2.w;
        s.t = s1.t + temp + s2.t;
        s.c = s1.c + s2.w * (s1.t + temp) + s2.c;
        s.first = s1.first; 
        s.last = s2.last;

        return s;
    }
} Subsequence;

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
void initiateNode(int i, Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    subseq_matrix[i][i].w = (i > 0);
    subseq_matrix[i][i].t = 0;
    subseq_matrix[i][i].c = 0;
    subseq_matrix[i][i].first = s->sequence[i];
    subseq_matrix[i][i].last = s->sequence[i];
}
void UpdateAllSubseq(Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    int n = s->sequence.size();

    for(int i = 0; i < n; i ++) initiateNode(i, s, subseq_matrix);

    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);

    for(int i = n-1; i >= 0; i--)
        for(int j = i-1; j >= 0; j--)
            subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]);
}
void VerifySubseq(Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    vector<vector<Subsequence>> subseq_matrixc(subseq_matrix);
    UpdateAllSubseq(s, subseq_matrixc);

    for(int i = 0; i < subseq_matrix.size(); i++)
    {
        for(int j = 0; j < subseq_matrix.size(); j++)
        {
            cout<<setfill(' ')<<setw(5)<<subseq_matrixc[i][j].c-subseq_matrix[i][j].c<<" ";
        }
        cout<<endl;
    }
}

Solution Construct()
{
    double alpha = 0.01f * (rand()%26);
    Solution s = Solution();

    s.sequence = {1};

    vector<int> cl;
    for(int i = 2; i <= SIZE; i++) cl.push_back(i);
    int r = 1;

    while(!cl.empty())
    {
        sort(cl.begin(), cl.end(),
        [r](int a, int b){return dt->getDistance(r, a) < dt->getDistance(r, b);});

        int limitOfBest = ceil(alpha * cl.size());
        int c = cl[limitOfBest>0?rand()%limitOfBest : 0];
        
        s.sequence.push_back(c);
        r = c;
        cl.erase(find(cl.begin(),cl.end(),c));
    }
    s.sequence.push_back(1);
    return s;
}

bool bestImprovementSwap(Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    double bestCost = s->cost;
    int best_i, best_j;
    int n = s->sequence.size();

    for(int i = 1; i < n - 1; i++)
    {
        for(int j = i+2; j < n - 1; j++)
        {
            Subsequence subs = subseq_matrix[0][i-1];
            subs = Subsequence::Concatenate(subs, subseq_matrix[j][j]);
            subs = Subsequence::Concatenate(subs, subseq_matrix[i+1][j-1]);
            subs = Subsequence::Concatenate(subs, subseq_matrix[i][i]);
            subs = Subsequence::Concatenate(subs, subseq_matrix[j+1][n-1]);

            if(subs.c < bestCost)
            {
                bestCost = subs.c;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    if(bestCost < s->cost) 
    {
        std::swap(s->sequence[best_i], s->sequence[best_j]);
        s->cost = bestCost;

        initiateNode(best_i, s, subseq_matrix);
        initiateNode(best_j, s, subseq_matrix);

        for(int i = 0; i <= best_i; i++)
            for(int j = best_i + (i == best_i); j < n; j++)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);
        for(int i = best_i+1; i <= best_j; i++)
            for(int j = best_j + (i == best_j); j < n; j++)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);

        for(int i = best_j-1; i >= best_i; i--)
            for(int j = best_i - (i == best_i); j >= 0; j--)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]);
        for(int i = n-1; i >= best_j; i--)
            for(int j = best_j - (i == best_j); j >= 0; j--)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]);

        return  true;
    }
    else return false;
}
bool bestImprovement2Opt(Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    double bestCost = s->cost;
    int best_i, best_j;
    int n = s->sequence.size();

    for(int i = 0; i < s->sequence.size()-1; i++)
    {
        for(int j = i+2; j < s->sequence.size()-1; j++)
        {
            Subsequence subs = subseq_matrix[0][i];
            subs = Subsequence::Concatenate(subs, subseq_matrix[j][i+1]);
            subs = Subsequence::Concatenate(subs, subseq_matrix[j+1][n-1]);
            
            if(subs.c < bestCost)
            {
                bestCost = subs.c;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(bestCost < s->cost) 
    {
        reverse(s->sequence.begin() + best_i + 1, s->sequence.begin() + best_j + 1);
        s->cost = bestCost;

        for(int i = best_i + 1; i <= best_j; i++) initiateNode(i, s, subseq_matrix);

        for(int i = 0; i <= best_j; i++)
            for(int j = (i>best_i)?(i+1):(best_i+1); j < n; j++)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);
        for(int i = n-1; i > best_i; i--)
            for(int j = (i<=best_j)?(i-1):(best_j); j >= 0; j--)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]);

        return true;
    }
    else return false;
}
bool bestImprovementOrOpt(Solution *s, vector<vector<Subsequence>> &subseq_matrix, int nVertex)//TODO
{
    double bestCost = s->cost;
    int best_block, best_edge;
    int n = s->sequence.size();

    for(int block_start_index = 1; block_start_index < n-1 - (nVertex-1); block_start_index++)
    {
        int block_end_index = block_start_index + nVertex-1;
        
        for(int j = 0; j < n-1 - (nVertex+2); j++)
        {
            Subsequence subs;
            int edge_start_index = j;

            if(j >= block_start_index-1)
            {
                edge_start_index += nVertex+2;

                subs = subseq_matrix[0][block_start_index-1];
                subs = Subsequence::Concatenate(subs,subseq_matrix[block_end_index+1][edge_start_index]);
                subs = Subsequence::Concatenate(subs,subseq_matrix[block_start_index][block_end_index]);
                subs = Subsequence::Concatenate(subs,subseq_matrix[edge_start_index+1][n-1]);
            }
            else
            {
                subs = subseq_matrix[0][edge_start_index];
                subs = Subsequence::Concatenate(subs,subseq_matrix[block_start_index][block_end_index]);
                subs = Subsequence::Concatenate(subs,subseq_matrix[edge_start_index+1][block_start_index]);
                subs = Subsequence::Concatenate(subs,subseq_matrix[block_end_index+1][n-1]);
            }

            if(subs.c < bestCost)
            {
                bestCost = subs.c;
                best_block = block_start_index;
                best_edge = edge_start_index;
            }
        }
    }

    if(bestCost < s->cost) 
    {
        int higher_bound, lower_bound;
        if(best_block<best_edge)
        {
            rotate(
                s->sequence.begin()+best_block, 
                s->sequence.begin()+best_block + nVertex, 
                s->sequence.begin()+best_edge + 1);
            higher_bound = best_edge;
            lower_bound = best_block;
        }
        else 
        {
            rotate(
                s->sequence.begin()+best_edge + 1, 
                s->sequence.begin()+best_block, 
                s->sequence.begin()+best_block + nVertex);
            higher_bound = best_block;
            lower_bound = best_edge;
        }
        s->cost = bestCost;
        
        for(int i = lower_bound; i <= higher_bound; i++) initiateNode(i, s, subseq_matrix);

        for(int i = 0; i <= higher_bound; i++)
            for(int j = (i>=lower_bound)?(i+1):(lower_bound); j < n; j++)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j-1], subseq_matrix[j][j]);
        for(int i = n-1; i >= lower_bound; i--)
            for(int j = (i<=higher_bound)?(i-1):(higher_bound); j >= 0; j--)
                subseq_matrix[i][j] = Subsequence::Concatenate(subseq_matrix[i][j+1], subseq_matrix[j][j]);

        return  true;
    }
    else return false;
}

void LocalSearch(Solution *s, vector<vector<Subsequence>> &subseq_matrix)
{
    std::vector<int> NL = {1, 2, 3, 4, 5};
    bool improved = false;

    while(!NL.empty())
    {
        int n = rand() % NL.size();

        switch (NL[n]) 
        {
        case 1:
            improved = bestImprovementSwap(s, subseq_matrix);
            break;
        case 2:
            improved = bestImprovement2Opt(s, subseq_matrix);
            break;
        case 3:
            improved = bestImprovementOrOpt(s,subseq_matrix,1);
            break;
        case 4:
            improved = bestImprovementOrOpt(s,subseq_matrix,2);
            break;
        case 5:
            improved = bestImprovementOrOpt(s,subseq_matrix,3);
            break;
        }

        if(improved) 
        {
            NL = {1,2,3,4,5};
        }
        else NL.erase(NL.begin() + n);
    }
} 
Solution Pertubation(const Solution &s, vector<vector<Subsequence>> &subseq_matrix)//TODO
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
        vector<vector<Subsequence>> subseq_matrix(s.sequence.size(), 
            vector<Subsequence>(s.sequence.size()));
        UpdateAllSubseq(&s, subseq_matrix);
        s.cost = subseq_matrix[0][s.sequence.size()-1].c;

        Solution best = s;

        int iterIls = 0;

        while(iterIls <= maxIterIls)
        {
            LocalSearch(&s, subseq_matrix);
            
            if(s.cost < best.cost)
            {
                best = s;
                iterIls = 0;
            }

            //s = Pertubation(best, subseq_matrix);
            
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

    cout << "Begining ILS" << endl;
    clock_t before = clock();
    Solution s = ILS(10, SIZE<100?SIZE:100);

    float duration = (clock()-before);
    s.print("ILS Solution");
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;
    
    return 0;
}