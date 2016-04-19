//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#include "cardinalityestimation.hpp"
#include "hyperloglog.hpp"

#include <iostream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <functional>
#include <chrono>

using namespace std;

const int minP = 4;
const int maxP = 22;
const int maxPplusQ = 64;

const string seedFileName = "../data/hll/seed.dat";
const string cardFileName = "../data/cardinalities.dat";
const string dataPathName = "../data/hll/";
const string intersectDataPathName = "../data/hll_joint/";
const string jointCardFileName = intersectDataPathName + "joint_cardinalities.dat";

template <typename T> void printToFile(const string& fileName, const T* data, size_t numSketches, bool append) {
    ofstream file(fileName, (append)?(ofstream::out | ofstream::app):ofstream::out);
    for (size_t s = 0; s < numSketches; ++s) {
        if (s>0) {
            file << " ";
        }
        file << data[s];
    }
    file << endl;
}

vector<int> readCounts(string s) {
    stringstream stream(s);
    int i;
    vector<int> v;
    while(stream >> i) {
        v.push_back(i);
    }
    v.shrink_to_fit();
    return v;
}

string getFileName(const string& seed, int p) {
    stringstream s;
    s << dataPathName;
    s << seed;
    s << "_" << setfill('0') << setw(2) << dec << p;
    s << ".dat";
    return s.str();
}



string getOutputFileName(const long cardA, long cardB, long cardX, int p, int q) {
    stringstream s;
    s << setfill('0') << setw(11) << cardA;
    s << "_";
    s << setfill('0') << setw(11) << cardB;
    s << "_";
    s << setfill('0') << setw(11) << cardX;
    s << "_";
    s << setfill('0') << setw(2) << p;
    s << "_";
    s << setfill('0') << setw(2) << q;
    s << ".dat";
    return s.str();
}




int main(int argc, char* argv[])
{
    mt19937_64 rng(chrono::system_clock::now().time_since_epoch().count());

    
    vector<string> seeds;
    {
        ifstream seedFile(seedFileName);
        
        string line;
        while (getline(seedFile, line))
        {
            seeds.push_back(line);
        }
    }
    
    cout << "seedSize "  << seeds.size() << endl;
    
    vector<long> cardinalities;
    {
        ifstream cardFile(cardFileName);
        
        string line;
        while (getline(cardFile, line))
        {
            cardinalities.push_back(stol(line));
        }
    }
    
    cout << "cardinalitySize "  << cardinalities.size() << endl;
        
    const unsigned char p = 20;
    const unsigned char q = 44;
    
    #pragma omp parallel
    while(true) {
        
        uint_fast64_t seed;
        #pragma omp critical
        {
            seed = rng();
        }
        
        mt19937_64 rng2(seed);
        std::uniform_int_distribution<int> uni(0,cardinalities.size()-1);
    
        const int cardIdxX = uni(rng2);
        const int cardIdxA = uni(rng2);
        const int cardIdxB = uni(rng2);
        
        long trueCardA = cardinalities[cardIdxA];
        long trueCardB = cardinalities[cardIdxB];
        long trueCardX = cardinalities[cardIdxX];

        string outputFileName = getOutputFileName(trueCardA, trueCardB, trueCardX, p, q);
        ofstream outputFile(intersectDataPathName + outputFileName);
        
        for(size_t seedCounter = 0; seedCounter+2 < seeds.size(); seedCounter+=3) {
            
            const size_t seedCounterA = seedCounter + 0;
            const size_t seedCounterB = seedCounter + 1;
            const size_t seedCounterX = seedCounter + 2;
            
            vector<int> cA(q+2);
            vector<int> cB(q+2);
            vector<int> cX(q+2);
            {
                ifstream dataFile(getFileName(seeds[seedCounterA], p));
                string line;
                for(size_t lineIdx = 0; lineIdx <= cardIdxA; ++lineIdx) {
                    getline(dataFile, line);
                }
                vector<int> tmpC = readCounts(line);
                for (int i = 0; i < tmpC.size(); ++i) {
                    cA[min(i, q+1)] += tmpC[i];
                }
            }
            {
                ifstream dataFile(getFileName(seeds[seedCounterB], p));
                string line;
                for(size_t lineIdx = 0; lineIdx <= cardIdxB; ++lineIdx) {
                    getline(dataFile, line);
                }
                vector<int> tmpC = readCounts(line);
                for (int i = 0; i < tmpC.size(); ++i) {
                    cB[min(i, q+1)] += tmpC[i];
                }
            }
            {
                ifstream dataFile(getFileName(seeds[seedCounterX], p));
                string line;
                for(size_t lineIdx = 0; lineIdx <= cardIdxX; ++lineIdx) {
                    getline(dataFile, line);
                }
                vector<int> tmpC = readCounts(line);
                for (int i = 0; i < tmpC.size(); ++i) {
                    cX[min(i, q+1)] += tmpC[i];
                }
            }
            
            HyperLogLog hllA = HyperLogLog::createFromCounts(cA, rng2);
            HyperLogLog hllB = HyperLogLog::createFromCounts(cB, rng2);
            HyperLogLog hllX = HyperLogLog::createFromCounts(cX, rng2);
            
            HyperLogLog hll1 = HyperLogLog::merge(hllA, hllX);
            HyperLogLog hll2 = HyperLogLog::merge(hllB, hllX);
            
            const vector<int> jointStatistic = HyperLogLog::getJointStatistic(hll1, hll2);
            
            assert(jointStatistic.size() == (q+2)*5);
            
            int pos = 0;
            for(int pos = 0; pos < (q+2)*5; ++pos) {
                outputFile << jointStatistic[pos] << " ";
            }
            outputFile << '\n';
            
        }
        #pragma omp critical
        {
            ofstream jointCardFile(jointCardFileName, std::ios::out | std::ios::app);
            jointCardFile << outputFileName << endl;
        }

    }
}
