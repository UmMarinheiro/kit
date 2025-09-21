#include "vectorUtils.h"
#include <stdexcept>
using namespace std;

vector<double> add(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing vector addition");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = a[i] + b[i];
            
    return c;
}
vector<double> subtr(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing vector subtraction");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = a[i] - b[i];
            
    return c;
}
vector<double> multiply(double a, const vector<double>& b)
{
    vector<double> c(b.size());
    for(int i = 0; i < b.size(); i++) c[i] = a * b[i];
            
    return c;
}
double dot(const vector<double>& a, const vector<double>& b)
{
    if(a.size() != b.size()) throw runtime_error("Size mismatch when performing dot product");

    double c = 0;
    for(int i = 0; i < a.size(); i++) c += a[i]*b[i];
            
    return c;
}
vector<double> apply(const vector<vector<double>>& a, vector<double>& b)
{
    if(a[0].size() != b.size()) throw runtime_error("Size mismatch when performing matrix-vector multiplication");

    vector<double> c(a.size());
    for(int i = 0; i < a.size(); i++) c[i] = dot(a[i], b);
            
    return c;
}

bool isLesserOrEqualTo0(const vector<double>& a)
{
    for(double i : a) if(i > 0) return false;
    return true;
}