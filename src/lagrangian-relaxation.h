#pragma once

#include <limits>
#include "Kruskal.h"
#include "vectorUtils.h"

vector<double> packCostsForLagrangean();
vector<vector<double>> getA();
vector<double> getOptimal(const vector<double>& lambda, const vector<pair<int,int>>& forbidden_arcs, double* vObj = nullptr);
vector<double> SolveLagrangianDual(double UB, const vector<double>& c,
  const vector<vector<double>>& a, const vector<double>& b,
  const vector<pair<int, int>>& forbidden_arcs);
