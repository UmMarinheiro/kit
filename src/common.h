#pragma once

#include <string>
#include <vector>

extern int n;
double getMatrixCost(int i, int j);
void initializeData(char* file);
void freeData();

typedef struct Solution
{
    std::vector<int> sequence;
    double cost;
    void print() const;
    void print(std::string name) const;
}Solution;

void updateCost(Solution& s);
