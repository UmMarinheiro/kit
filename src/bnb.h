#include "common.h"
#include <csignal>
#include <list>
#include <queue>

template<typename T>
class BranchingMethod
{
public:
  virtual void push(T toPush) = 0;
  virtual T pull() = 0;
  virtual bool empty() = 0;
protected:
  BranchingMethod(){}
};

template<typename T>
class DFS : public BranchingMethod<T>
{
  std::list<T> tree;
public:
  void push(T toPush)
  {
    tree.push_back(toPush);
  }
  T pull()
  {
    T ans = tree.back();
    tree.pop_back();
    return ans;
  }
  bool empty() {return tree.empty();}
};

template<typename T>
class BFS : public BranchingMethod<T>
{
  std::list<T> tree;
public:
  void push(T toPush)
  {
    tree.push_back(toPush);
  }
  T pull()
  {
    T ans = tree.front();
    tree.pop_front();
    return ans;
  }
  bool empty() {return tree.empty();}
};

template<typename T>
class Priority : public BranchingMethod<T>
{
  std::priority_queue<T> tree;
public:
  void push(T toPush)
  {
    tree.push(toPush);
  }
  T pull()
  {
    T ans = tree.top();
    tree.pop();
    return ans;
  }
  bool empty(){return tree.empty();}
};

typedef struct Node
{
  double cost = 0;
  std::vector<double> lambda;
  std::vector<double> sol;

  std::vector<std::pair<int, int>> forbidden_arcs;

  bool operator<(Node other) const{return cost > other.cost;}
  Node(double UB, const std::vector<double> &c, const std::vector<std::vector<double>> &a, const std::vector<double> &b, const std::vector<std::pair<int, int>> &fa);
}Node;

extern BranchingMethod<Node>* treePtr;

Solution bnb();
