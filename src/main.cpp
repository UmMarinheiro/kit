#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include "data.h"
#include "Kruskal.h"
#include "vectorUtils.h"
using namespace std;

#define EMIN 1e-5
#define KMAX 30

vector<double> getOptimal();

vector<double> SolveLagrangianDual(double UB, const vector<double>& c, const vector<vector<double>>& a, vector<double>& b)
{
    vector<double> bestlambda(c.size(), 0);
    vector<double> lambda(c.size(), 0);
    
    double bestw = numeric_limits<double>().lowest();
    
    vector<double> currlambda(c.size(), 0);
    double e = 1;
    double k = 0;
    
    vector<double> b_Ax;
    double w;
    do {
        vector<double> x = getOptimal();

        b_Ax = subtr(b, apply(a, x));
        w = dot(c, x) + dot(lambda, b_Ax);

        if(w > bestw)
        {
            bestw = w;
            bestlambda = lambda;
            k = 0;
        }
        else 
        {
            k++;
            if(k>KMAX)
            {
                k = 0;
                e /= 2;
            }    
        }
        double u = e*(UB - w)/dot(b_Ax, b_Ax);
        lambda = add(lambda,  multiply(u, b_Ax));
    }while (!(e < EMIN || (isLesserOrEqualTo0(b_Ax) && dot(lambda, b_Ax) == 0) || w >= UB)); 

    return bestlambda;
}

int main(int argc, char** argv) 
{
	initializeData(argv[1]);
	initializeCost();

	cout << "Running TSP..." << endl;
    clock_t before = clock();

    float duration = (clock()-before);
    s.print("TSP Solution");
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;
	
	freeData();

	return 0;
}
