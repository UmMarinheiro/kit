#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include "data.h"
#include "Kruskal.h"
using namespace std;

#define EMIN 1e-5
#define KMAX 30

vector<double> add(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing vector addition");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = a[i] + b[i];
            
    return c;
}
vector<double> subtr(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing vector subtraction");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = a[i] - b[i];
            
    return c;
}
vector<double> multiply(double a, const vector<double>& b)
{
    vector<double> c(b.size());
    for(int i = 0; i < b.size(); i++) c[i] = a * b[i];
            
    return c;
}
double dot(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing dot product");

    double c = 0;
    for(int i = 0; i < a.size(); i++) c += a[i]*b[i];
            
    return c;
}
vector<double> apply(const vector<vector<double>>& a, vector<double>& b)
{
    if(a[0].size() != b.size()) throw runtime_error("Size mismatch when performing matrix-vector multiplication");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = dot(a[i], b);
            
    return c;
}

bool isLesserOrEqualTo0(const vector<double>& a)
{
    for(double i : a) if(i > 0) return false;
    return true;
}
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
