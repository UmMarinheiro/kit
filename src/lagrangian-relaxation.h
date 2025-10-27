#include <limits>
#include "Kruskal.h"
#include "vectorUtils.h"

vector<double> packCostsForLagrangean();
vector<vector<double>> getA();
vector<double> SolveLagrangianDual(double UB, const vector<double>& c,
  const vector<vector<double>>& a, const vector<double>& b);
