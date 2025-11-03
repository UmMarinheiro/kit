#include "lagrangian-relaxation.h"
#include "common.h"
#include "defines.h"

#define INFINITE 9999999

// Inclui todos os nós menos o [0][0]
// e  adiciona o custo dos penalizadores
vector<vector<double>> packCostsForKruskal(const vector<double>& lambda, const vector<pair<int,int>>& forbidden_arcs)
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
    for(const auto &arc : forbidden_arcs)
    {
        if(arc.first == 0 || arc.second == 0) continue;

        c[arc.first-1][arc.second-1] = INFINITE;
    }
    return c;
}
// Transforma a matrix de custos n x n em um vetor n2 x 1
vector<double> packCostsForLagrangean()
{
    vector<double> c(n*n, 0);
    for(int i = 0; i < n; i++)
        for(int j = i+1; j < n; j++)
        {
            c[i*n + j] = getMatrixCost(i, j);
        }
    return c;
}
bool isForbiden(pair<int,int> arc, vector<pair<int, int>> fa)
{
    for(auto &current : fa)
        if(current.first == arc.first && current.second == arc.second) return true;
    return false;
}
vector<double> getOptimal(const vector<double>& lambda, const vector<pair<int,int>>& forbidden_arcs, double* vObj)
{
    vector<vector<double>> costs = packCostsForKruskal(lambda, forbidden_arcs);
    Kruskal kruskal(costs);
    if(vObj != nullptr) *vObj = kruskal.MST(n-1);
    else kruskal.MST(n-1);
    
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
        if(isForbiden({0,i}, forbidden_arcs)) continue;

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
    x[node1*n + 0] = 1;
    x[0*n + node2] = 1;
    x[node2*n + 0] = 1;

    if(vObj != nullptr) *vObj += lowestDistToZero + secondLowestDistToZero;

    return x;
}

vector<double> SolveLagrangianDual(double UB, const vector<double>& c, const vector<vector<double>>& a, const vector<double>& b,
    const vector<pair<int, int>>& forbidden_arcs, vector<double> lambda)
{
    vector<double> bestlambda(n, 0);
    
    double bestw = numeric_limits<double>().lowest();
    
    double e = 1;
    double k = 0;

    while(true)
    {
        double w;

        vector<double> x = getOptimal(lambda, forbidden_arcs, &w);
        vector<double> b_Ax = subtr(b, apply(a, x));
        
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
        
        lambda = add(lambda, multiply(u, b_Ax));

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
