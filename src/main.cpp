#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>
#include "data.h"
#include "Kruskal.h"
#include "vectorUtils.h"
using namespace std;

#define EMIN 1e-5
#define KMAX 30

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

// Inclui todos os nós menos o [0][0]
// e  adiciona o custo dos penalizadores
vector<vector<double>> packCostsForKruskal(const vector<double>& lambda)
{
    vector<vector<double>> c(n-1, vector<double>(n-1));
    for(int i = 0; i < n-1; i++)
    {
        for(int j = 0; j < n-1; j++)
        {
            if(i < j) c[i][j] = getMatrixCost(i+1, j+1);
            else c[i][j] = getMatrixCost(j+1, i+1);
            c[i][j] -= lambda[i+1] + lambda[j+1];
        }
    }
    return c;
}
// Transforma a matrix de custos n x n em um vetor n2 x 1
vector<double> packCostsForLagrangean()
{
    vector<double> c(n*n);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
        {
            c[i*n + j] = getMatrixCost(i, j);
        }
    return c;
}

vector<double> getOptimal(const vector<double>& lambda)
{
    vector<vector<double>> costs = packCostsForKruskal(lambda);
    Kruskal kruskal(costs);
    kruskal.MST(n-1);
    
    vector<double> x(n*n, 0);
    auto edges = kruskal.getEdges();

    for(auto [i,j] : edges)
    {
         x[(i+1)*n + (j+1)] = 1;  
         x[(j+1)*n + (i+1)] = 1;  
    } 

    double lowestDistToZero = numeric_limits<double>().max();
    int node1 = 0;
    double secondLowestDistToZero = numeric_limits<double>().max();
    int node2 = 0;

    for(int i = 1; i < n; i++) 
    {
        double currentCost = getMatrixCost(0, i) - lambda[0] - lambda[i];
        if(currentCost < lowestDistToZero) 
        {
            node2 = node1;
            secondLowestDistToZero = lowestDistToZero;
            node1 = i;
            lowestDistToZero = currentCost;
        }
        else if(currentCost < secondLowestDistToZero)
        {
            node2 = i;
            secondLowestDistToZero = currentCost;
        }
    }
    x[0*n + node1] = 1;
    x[node1*n] = 1;
    x[0*n + node2] = 1;
    x[node2*n + 0] = 1;

    return x;
}

vector<double> SolveLagrangianDual(double UB, const vector<double>& c, const vector<vector<double>>& a, vector<double>& b)
{
    vector<double> bestlambda(n, 0);
    vector<double> lambda(n, 0);
    
    double bestw = numeric_limits<double>().lowest();
    
    double e = 1;
    double k = 0;

    while(true)
    {
        vector<double> x = getOptimal(lambda);
        vector<double> b_Ax = subtr(b, apply(a, x));
        
        double w = dot(c, x) + dot(lambda, b_Ax);
        
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

        if(e < EMIN) break;
        if(isLesserOrEqualTo0(b_Ax) && dot(lambda, b_Ax)) break;
        if(w >= UB) break;
    }

    return bestlambda;
}

vector<vector<double>> getA()
{
    vector<vector<double>> a(n, vector<double>(n*n, 0));

    for(int i = 0; i < n; i++)
    {
        for(int j = i+1; j < n; j++) a[i][i*n + j] = 1;
        for(int j = i-1; j >= 0; j--) a[i][j*n + i] = 1;
    }

    return a;
}

int main(int argc, char** argv) 
{
    initializeData(argv[1]);

    cout << "Running TSP..." << endl;
    clock_t before = clock();

    vector<double> c = packCostsForLagrangean();
    vector<vector<double>> a = getA();
    vector<double> b(n, 2);
    auto dd = SolveLagrangianDual(148, c, a, b);
    cout<<endl;
    for(int i = 0; i < n; i++) cout<<dd[i]<<" ";
    cout<<std::endl;
                        
    float duration = (clock()-before);
    cout << "Took " << (float)duration/CLOCKS_PER_SEC << " seconds" << endl;

    freeData();

    return 0;
}
