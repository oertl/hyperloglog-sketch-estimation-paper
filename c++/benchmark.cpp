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
#include <memory>
#include <chrono>

using namespace std;

const int minP = 4;
const int maxP = 22;
const int maxPplusQ = 64;


const int numLoops = 100;

const int numData = 1000;

const string seedFileName = "../data/hll/seed.dat";
const string cardFileName = "../data/cardinalities.dat";
const string dataPathName = "../data/hll/";

const int p1 = 12;
const int q1 = 52;
const string resultsFileName1 = "../data/benchmark_results_12_52.dat";
const int p2 = 12;
const int q2 = 20;
const string resultsFileName2 = "../data/benchmark_results_12_20.dat";

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

void run(const int p, const int q, const string& resultsFileName) {
    
    
    mt19937_64 rng(0);

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
    
    assert(numData <= seeds.size());
    
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
    
    
    std::vector<unique_ptr<ifstream>> files;
    files.reserve(numData);
    for (size_t seedCounter = 0; seedCounter < numData; ++seedCounter) {
        files.emplace_back(unique_ptr<ifstream>(new ifstream(getFileName(seeds[seedCounter], p))));
    }
    
    size_t dataSize = numLoops * numData;

    double sum = 0.;
    
    ofstream outputFile(resultsFileName);

    for(long cardinality : cardinalities) {

        // load data into memory
        vector<vector<int>> counts(dataSize);
        for(size_t seedCounter = 0; seedCounter < numData; ++seedCounter) {
        
            ifstream& dataFile = (*files[seedCounter]);
            
            string line;
            getline(dataFile, line);
            
            vector<int> tmpC = readCounts(line);
            vector<int> c(q+2);
            for (int i = 0; i < tmpC.size(); ++i) {
                c[min(i, q+1)] += tmpC[i];
            }
            for (int i = 0; i < numLoops; ++i) {
                counts[seedCounter + i *numData] = c;
            }
        }
        
        /*vector<HyperLogLog> hyperLogLogs;
        hyperLogLogs.reserve(dataSize);
        for (auto c : counts) {
            hyperLogLogs.emplace_back(HyperLogLog::createFromCounts(c, rng));
        }*/
        
        double maxLikelihoodWithRegisterScan = 0.;
        double maxLikelihoodWithoutRegisterScan = 0.;
        double flajoletWithRegisterScan = 0.;
        double flajoletWithoutRegisterScan = 0.;
        double flajoletCorrectedWithoutRegisterScan = 0.;
        long flajoletCorrectedSmallCorrectionIterationsSum = 0;
        long flajoletCorrectedLargeCorrectionIterationsSum = 0;
        long outerLoopIterationsCountSum = 0;
        long innerLoop1IterationsCountSum = 0;
        long innerLoop2IterationsCountSum = 0;
        long logEvaluationCountSum = 0;

        /*{
            auto start = std::chrono::system_clock::now();
            for(auto h : hyperLogLogs) {
                sum += maxLikelihoodEstimate(h.getCounts());
                
            }
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start);
            maxLikelihoodWithRegisterScan = elapsed.count()/(double)dataSize;
        }*/
        
        {
            int outerLoopIterationsCount;
            int innerLoop1IterationsCount;
            int innerLoop2IterationsCount;
            int logEvaluationCount;
            int kMin;
            int kMax;
            
            
            auto start = std::chrono::system_clock::now();
            for(auto c : counts) {
                sum += maxLikelihoodEstimate(c, outerLoopIterationsCount, innerLoop1IterationsCount, innerLoop2IterationsCount, logEvaluationCount, kMin, kMax);
                outerLoopIterationsCountSum += outerLoopIterationsCount;
                innerLoop1IterationsCountSum += innerLoop1IterationsCount;
                innerLoop2IterationsCountSum += innerLoop2IterationsCount;
                logEvaluationCountSum += logEvaluationCount;
            }
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start);
            maxLikelihoodWithoutRegisterScan = elapsed.count()/(double)dataSize;
        }
        
        {
            int smallCorrectionIterations;
            int largeCorrectionIterations;
            
            auto start = std::chrono::system_clock::now();
            for(auto c : counts) {
                sum += flajoletRawEstimateCorrected(c, smallCorrectionIterations, largeCorrectionIterations);
                flajoletCorrectedSmallCorrectionIterationsSum += smallCorrectionIterations;
                flajoletCorrectedLargeCorrectionIterationsSum += largeCorrectionIterations;
            }
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start);
            flajoletCorrectedWithoutRegisterScan = elapsed.count()/(double)dataSize;
        }
            
        /*{
            auto start = std::chrono::system_clock::now();
            for(auto h : hyperLogLogs) {
                sum += flajoletEstimate(h.getCounts());
                
            }
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start);
            flajoletWithRegisterScan = elapsed.count()/(double)dataSize;
        }*/
        
        /*{
            auto start = std::chrono::system_clock::now();
            for(auto c : counts) {
                sum += flajoletEstimate(c);
            }
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - start);
            flajoletWithoutRegisterScan = elapsed.count()/(double)dataSize;
        }*/
        
        outputFile << cardinality;
        //outputFile << " " << maxLikelihoodWithRegisterScan;
        outputFile << " " << maxLikelihoodWithoutRegisterScan;
        //outputFile << " " << flajoletWithRegisterScan;
        //outputFile << " " << flajoletWithoutRegisterScan;
        outputFile << " " << outerLoopIterationsCountSum/(double)dataSize;
        outputFile << " " << innerLoop1IterationsCountSum/(double)dataSize;
        outputFile << " " << innerLoop2IterationsCountSum/(double)dataSize;
        outputFile << " " << logEvaluationCountSum/(double)dataSize;
        outputFile << " " << flajoletCorrectedWithoutRegisterScan;
        outputFile << " " << flajoletCorrectedSmallCorrectionIterationsSum/(double)dataSize;
        outputFile << " " << flajoletCorrectedLargeCorrectionIterationsSum/(double)dataSize;
        outputFile << endl;
    }
    
    cout << sum << endl; // in order to avoid compiler optimizations
}

int main(int argc, char* argv[])
{
    run(p1, q1, resultsFileName1);
    run(p2, q2, resultsFileName2);
    
    return 0;
}

