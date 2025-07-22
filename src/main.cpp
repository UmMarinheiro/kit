#include <cstddef>
#include <iostream>

#include <ctime>
#include <vector>
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
    Solution s = {{1,2,3,4,5,6,1}, 0};//ILS(50, SIZE/(1 + (SIZE>=150)));
    updateCost(s);
    float duration = (clock()-before);
    s.print("ILS Solution");
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}