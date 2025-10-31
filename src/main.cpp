#include <cmath>
#include <cstdio>
#include <iostream>
#include "common.h"
#include "bnb.h"

using namespace std;

int main(int argc, char** argv) 
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
