#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include "lagrangian-relaxation.h"
#include "common.h"
#include "grasp.h"

using namespace std;

int main(int argc, char** argv) 
{
    initializeData(argv[1]);

    cout << "Running TSP..." << endl;
    clock_t before = clock();

    vector<double> c = packCostsForLagrangean();
    vector<vector<double>> a = getA();
    vector<double> b(n, 2);
    vector<pair<int, int>> forbidden_arcs = {};
    double vObj;
    
    Solution s = Construct();

    cout<<"UB: "<< s.cost<<endl; 

    auto dd = SolveLagrangianDual(s.cost, c, a, b, forbidden_arcs);
    getOptimal(dd, forbidden_arcs, &vObj);

    cout<<endl<<endl;
    cout<<"vObj: " << vObj<<endl;
    cout<<endl;
    for(int i = 0; i < n; i++) cout<<dd[i]<<" ";
    cout<<std::endl;
                        
    float duration = (clock()-before);
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;

    freeData();

    return 0;
}
