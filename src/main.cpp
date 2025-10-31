#include <cmath>
#include <cstdio>
#include <iostream>
#include "common.h"
#include "bnb.h"

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

// Inclui todos os nós menos o [0][0]
// e  adiciona o custo dos penalizadores
vector<vector<double>> packCostsForKruskal(const vector<double>& lambda)
{
    initializeData(argv[1]);

    if(argc < 3)
    {
        treePtr = new DFS<Node>(); 
        cout<<"Mode set to DFS"<<endl;
    }
    else
    {
        if((string)argv[2] == "DFS")
        {
            treePtr = new DFS<Node>();
            cout<<"Mode set to DFS"<<endl;
        }
        else if((string)argv[2] == "BFS")
        {
            treePtr = new BFS<Node>();
            cout<<"Mode set to BFS"<<endl;
        }
        else if((string)argv[2] == "Lowest")
        {
            treePtr = new Priority<Node>();
            cout<<"Mode set to Lowest"<<endl;
        }
        else
        {
            treePtr = new DFS<Node>();
            cout<<"Mode set to DFS"<<endl;
        }
    }

    cout << "Running TSP..." << endl;
    clock_t before = clock();

    Solution sol = bnb();

    sol.print();
                        
    float duration = (clock()-before);
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;

    freeData();

    return 0;
}
