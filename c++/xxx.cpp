#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <cstdlib>

using namespace std;


double sigma(int c, int m, int& numIterations) {
        numIterations = 0;
        if (c==m) return std::numeric_limits<double>::infinity();
        double powerOfX = static_cast<double>(c)/static_cast<double>(m);
        double result = powerOfX;
        double previousResult;
        double powerOf2 = 1;
        do {
            numIterations += 1;
            powerOfX *= powerOfX;
            previousResult = result;
            result += powerOfX * powerOf2;
            powerOf2 += powerOf2;
        } while(previousResult < result);
        return m*result;
    }
double tau(int c, int m, int& numIterations) {
	numIterations = 0;
	double result = static_cast<double>(c)/static_cast<double>(m);
	double powerOfX  = 1 - result;
	double previousResult;
	double powerOf2 = 0.5;
	do {
		numIterations += 1;
		powerOfX = std::sqrt(powerOfX);
		previousResult = result;
		result -= (1 - powerOfX)*powerOf2;
		powerOf2 *= 0.5;
	} while(previousResult > result);
	return m*result;
}

double tau2(int c, int m, int& numIterations) {
	numIterations = 0;
	
	double powerOfX = static_cast<double>(m-c)/static_cast<double>(m);
	
	double result = 0;
	double previousResult;
	double powerOf2 = 1.0;
	do {
		numIterations += 1;
		powerOfX = std::sqrt(powerOfX);
		previousResult = result;
		result += (1 - powerOfX)*powerOfX*powerOf2;
		powerOf2 *= 0.5;
	} while(previousResult < result);
	return m*result;
}

int main(int argc, char* argv[]) {
	int m = 4096;

	for (int c =0; c <= m; ++c) {
		
		int numIterations;
		int numIterations2;
		
		double t = sigma(c, m, numIterations);
		double t2 = tau2(c, m, numIterations2);
		
		cout << c << " " << t << " " << t2 << " " << numIterations << " " << numIterations2 << endl;
	}


	return 0;
}
