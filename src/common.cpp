#include "data.h"

Data * dt;
int n;

double getMatrixCost(int i, int j)
{
    // static double distCost[5][5] =
    // {
    //     {0, 30, 26, 50, 40},
    //     {0,  0, 24, 40, 50},
    //     {0,  0,  0, 24, 26},
    //     {0,  0,  0,  0, 30},
    //     {0,  0,  0,  0,  0}
    // };
    // return distCost[i][j];
    return dt->getDistance(i, j);
}

void initializeData(char* file)
{
    cout << "Reading " << file << " ..." << endl;
    dt = new Data(2, file);
    dt->readData();
    n = dt->getDimension();
    // n = 5;
    cout << "Succesfully read " << file << " !" << endl;
}
void freeData() {delete dt;}


