//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#include "hyperloglog.hpp"

#include <fstream>
#include <random>
#include <iostream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <array>
#include <memory>

using namespace std;

const int maxP = 22;
const int minP = 4;

const string outputFolder = "../data/hll/";
const string seedFileName = "../data/hll/seed.dat";
const string cardFileName = "../data/cardinalities.dat";

string getFileName(const uint_fast64_t seed, int q) {
	stringstream s;
	s << outputFolder;
	s << hex <<  setfill('0') << setw(16) << seed;
	s << "_" << setw(2) << dec << q;
	s << ".dat";
	return s.str();
}

int main(int argc, char* argv[])
{
	// read cardinalities
	ifstream cardFile(cardFileName);
	vector<long> cardinalities;
	string line;
	while (getline(cardFile, line))
	{
		cardinalities.push_back(stol(line));
	}
	
	if (cardinalities.empty()) {
		return 1;
	}
	
	mt19937_64 rng(chrono::system_clock::now().time_since_epoch().count());
	
	#pragma omp parallel
	while(true) {
		
		uint_fast64_t seed;
		#pragma omp critical
		{
			seed = rng();
		}
		
		mt19937_64 rng2(seed);
		
		long lastCardinality = 0;
	
		HyperLogLog hll(maxP, 64-maxP);
		
		vector<vector<vector<int>>> data(maxP-minP+1);
		
		for(size_t cardinalityCounter = 0; cardinalityCounter < cardinalities.size(); cardinalityCounter+=1) {
			long cardinality = cardinalities[cardinalityCounter];
			for (long e = lastCardinality; e < cardinality; ++e) {
				hll.add(rng2());
			}
			HyperLogLog hll2 = hll;
			data[maxP-minP].push_back(hll2.getCounts());
			for (int p = maxP-1; p >= minP; --p) {
				hll2 = hll2.reduce(p, 64-p);
				data[p-minP].push_back(hll2.getCounts());
			}
			
			lastCardinality = cardinality;
		}
		
		for (int p = minP; p <= maxP; ++p) {
			ofstream f(getFileName(seed, p));
			for (const vector<int>& counts : data[p-minP]) {
				bool b = false;
				for (int count : counts) {
					if (b) {
						f << " ";
					}
					b = true;
					f << count;
				}
				f << '\n';
			}
		}
		
		#pragma omp critical
		{
			ofstream seedFile(seedFileName, std::ios::out | std::ios::app);
			seedFile << hex << setfill('0') << setw(16)  << seed << endl;
		}
	}
}
