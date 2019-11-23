//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace std;

const double factor = 1.01;
const long maxCardinality = 50000000000;

int main(int argc, char* argv[])
{
    vector<long> cardinalities;
    
    long cardinality = maxCardinality;
    while(cardinality >= 0) {
        cardinalities.push_back(cardinality);
        cardinality = min(lround(cardinality/factor), cardinality-1);
    }
    
    sort(cardinalities.begin(), cardinalities.end());
    auto last = std::unique(cardinalities.begin(), cardinalities.end());
    cardinalities.erase(last, cardinalities.end()); 
    
    for(long c : cardinalities) {
        cout << c << endl;
    }
}
