#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include "common.h"
#include "grasp.h"
#include "bnb.h"

using namespace std;

int main(int argc, char** argv) 
{
    initializeData(argv[1]);

    cout << "Running TSP..." << endl;
    clock_t before = clock();

    Solution sol = bnb();

    sol.print();
                        
    float duration = (clock()-before);
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;

    freeData();

    return 0;
}
