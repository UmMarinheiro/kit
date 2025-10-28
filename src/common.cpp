#include "data.h"
#include "common.h"
#include "defines.h"
#include <vector>

Data * dt;
int n;

double getMatrixCost(int i, int j) 
{
    #ifdef TESTCASE
    static double distCost[5][5] =
    {
        {0, 30, 26, 50, 40},
        {0,  0, 24, 40, 50},
        {0,  0,  0, 24, 26},
        {0,  0,  0,  0, 30},
        {0,  0,  0,  0,  0}
    };
    return distCost[i][j];
    #else
    return dt->getDistance(i, j);
    #endif
}

void initializeData(char* file)
{
    cout << "Reading " << file << " ..." << endl;
    dt = new Data(2, file);
    dt->readData();
    #ifdef TESTCASE
    n = 5;
    #else
    n = dt->getDimension();
    #endif
    cout << "Succesfully read " << file << " !" << endl;
}
void freeData() {delete dt;}

void Solution::print() const
{
    for(int i = 0; i < sequence.size() - 1; i++)
        cout << sequence.at(i) << " -> ";
    cout << sequence.back() << std::endl;

    cout << "Cost: " << cost << std::endl;
}
void Solution::print(string name) const {cout<<name<<"="<<endl; print(); }

void updateCost(Solution &toUpdate)
{
    toUpdate.cost = 0;
    for(int i = 0; i < toUpdate.sequence.size()-1; i++) 
        toUpdate.cost += dt->getDistance(toUpdate.sequence[i], toUpdate.sequence[i+1]);
}
