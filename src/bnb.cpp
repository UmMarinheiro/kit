#include "common.h"
#include "grasp.h"
#include "lagrangian-relaxation.h"
#include <cstdio>
#include "bnb.h"

BranchingMethod<Node> *treePtr; 

Node::Node(double UB, const vector<double>& c, const vector<vector<double>>& a, const vector<double>& b, const vector<pair<int,int>>& fa)
{
  forbidden_arcs = fa;
  lambda = SolveLagrangianDual(UB, c, a, b, forbidden_arcs);
  sol = getOptimal(lambda, forbidden_arcs);
  for(int i = 0; i < n; i++)
    for(int j = i+1; j < n; j++)
      cost += sol[i*n + j]*getMatrixCost(i, j);
}

Solution bnb()
{
  double ub = Construct().cost;
  vector<double> bestSol;

  const vector<double> c = packCostsForLagrangean();
  const vector<vector<double>> a = getA();
  const vector<double> b(n, 2);

  if(treePtr == nullptr)
  {
    std::cout << "Branching method not defined!" << endl;
    throw; 
  }

  treePtr->push(Node(ub, c, a, b, {}));
  while(!treePtr->empty())
  {
    Node current = treePtr->pull();

    if(current.cost > ub) continue;

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

        treePtr->push(Node(ub, c, a, b, forbidden));
      }
    }
  }

  vector<int> sequence = {0};
  while(sequence.size() != n+1)
  {
    for(int i = 0; i < n; i++)
    {
      if(bestSol[sequence.back()*n + i] == 0) continue;
      
      bestSol[sequence.back()*n + i] = 0;
      bestSol[i*n + sequence.back()] = 0;
      sequence.push_back(i);
      break;
    }
  }

  return (Solution){sequence,ub};
}
