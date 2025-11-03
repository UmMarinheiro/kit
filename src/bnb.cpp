#include "common.h"
#include "grasp.h"
#include "lagrangian-relaxation.h"
#include <cstdio>
#include <vector>
#include "bnb.h"

BranchingMethod<Node> *treePtr; 

Node::Node(double UB, const vector<double>& c, const vector<vector<double>>& a, const vector<double>& b,
           const vector<pair<int,int>>& fa, const vector<double>& startingLambda)
{
  forbidden_arcs = fa;
  lambda = SolveLagrangianDual(UB, c, a, b, forbidden_arcs, startingLambda);
  sol = getOptimal(lambda, forbidden_arcs, &cost);
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

  treePtr->push(Node(ub, c, a, b, {}, vector<double>(n, 0)));
  while(!treePtr->empty())
  {
    Node current = treePtr->pull();

    cout << "Depth:" << current.forbidden_arcs.size()/2 << endl
      << "Cost: " << current.cost << "/" << ub << endl; 

    if(current.cost > ub - 1e-5) continue;

    int bestIndex = 0;
    int bestDegree = 0;
    bool feasible = true;

    for(int i = 0; i < n; i++)
    {
      int degree = 0;
      for(int j = 0; j < n; j++) degree += current.sol[i*n + j];
      if(degree > bestDegree)
      {
        bestIndex = i;
        bestDegree = degree;
      }
      if(degree != 2) feasible = false;
    }

    if(feasible)
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

        treePtr->push(Node(ub, c, a, b, forbidden, current.lambda));
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
