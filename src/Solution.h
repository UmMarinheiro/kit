#pragma once 

#include <vector>

typedef struct Solution
{
    std::vector<int> sequence;
    double cost;

    Solution(std::vector<int> s = {},double c = 0.0);

    const void print();
    const void print(char* name);
}Solution;
