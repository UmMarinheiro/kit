#pragma once
#include <vector>

std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> subtr(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> multiply(double a, const std::vector<double>& b);
double dot(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> apply(const std::vector<std::vector<double>>& a, std::vector<double>& b);
bool isLesserOrEqualTo0(const std::vector<double>& a);