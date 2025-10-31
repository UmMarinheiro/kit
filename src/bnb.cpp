#include "common.h"
#include "grasp.h"
#include "lagrangian-relaxation.h"
#include <cstdio>
#include <list>

typedef struct Node
{
  double cost = 0;
  vector<double> sol;
  vector<double> lambda;
  
  vector<pair<int, int>> forbidden_arcs;

  Node(double UB, const vector<double> &c, const vector<vector<double>> &a, const vector<double> &b, const vector<pair<int, int>> &fa)
  {
    this->forbidden_arcs = fa;
    lambda = SolveLagrangianDual(UB, c, a, b, fa);
    sol = getOptimal(lambda, fa);

    for(int i = 0; i < n; i++)
      for(int j = i+1; j < n; j++)
        cost += sol[i*n+j]*getMatrixCost(i, j);
  }
}node;

Solution bnb()
{
  double ub = Construct().cost;
  vector<double> bestSol;

  const vector<double> c = packCostsForLagrangean();
  const vector<vector<double>> a = getA();
  const vector<double> b(n, 2);

  list<Node> nodes = {Node(ub, c, a, b, {})};
  while(!nodes.empty())
  {
    auto itt = nodes.begin();
    node &current = *itt;

    cout<<"ub: "<< ub << " ccost: " << current.cost << endl;
    for(auto &i : current.forbidden_arcs) cout<<"[" << i.first << ", " << i.second << "], ";
    cout<<endl;

    if(current.cost > ub)
    {
      nodes.erase(itt);
      continue;
    }

    int bestIndex = 0;
    int bestDegree = 0;

    for(int i = 0; i < n; i++)
    {
      int degree = 0;
      for(int j = 0; j < n; j++) degree += current.sol[i*n + j];
      if(degree > bestDegree)
      {
        bestIndex = i;
        bestDegree = degree;
      }
    }

    if(bestDegree == 2)
    {
      ub = current.cost;
      bestSol = current.sol;
    }
    else
    {
      for(int i = 0; i < n; i++)
      {
        if(current.sol[bestIndex*n + i] == 0) continue;

        vector<pair<int, int>> forbidden = current.forbidden_arcs;

        forbidden.push_back({bestIndex, i});
        forbidden.push_back({i, bestIndex});

        nodes.push_back(Node(ub, c, a, b, forbidden));
      }
    }
    
    nodes.erase(itt);
  }

  vector<int> sequence = {0};
  while(sequence.size() != n+1)
  {
    for(int i = 0; i < n; i++)
    {
      if(bestSol[sequence.back()*n + i] == 1)
      {
        bestSol[sequence.back()*n + i] = 0;
        bestSol[i*n + sequence.back()] = 0;
        sequence.push_back(i);
        break;
      }
    }
  }

  return (Solution){sequence,ub};
}
