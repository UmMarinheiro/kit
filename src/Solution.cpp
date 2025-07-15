#include <vector>
#include <iostream>

#include "Solution.h"

Solution::Solution(std::vector<int> s,double c) : sequence(s), cost(c){}

const void Solution::print()
{
    for(int i = 0; i < sequence.size() - 1; i++)
    {
        if(sequence.at(i)<100)
        {
            std::cout << " ";
            if(sequence.at(i)<10) std::cout << " ";
        }
        
        std::cout << sequence.at(i) << " -> ";
    }
    std::cout << sequence.back() << std::endl;

    std::cout << "Cost: " << cost << std::endl;
}
const void Solution::print(char* name)
{
    std::cout << name << " = " << std::endl;
    print();
}
