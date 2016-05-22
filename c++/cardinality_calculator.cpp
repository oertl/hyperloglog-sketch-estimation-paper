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

using namespace std;

const int minP = 4;
const int maxP = 22;
const int maxPplusQ = 64;

const string seedFileName = "../data/hll/seed.dat";
const string cardFileName = "../data/cardinalities.dat";
const string dataPathName = "../data/hll/";

template <typename T> void printToFile(const string& fileName, const T* data, size_t numSketches, size_t numCardinalities) {
    ofstream file(fileName);
    for (size_t c = 0; c < numCardinalities; ++c) {
        for (size_t s = 0; s < numSketches; ++s) {
            if (s>0) {
                file << " ";
            }
            file << data[c*numSketches+s];
        }
        file << '\n';
    }
}

class parameter {
    int p;
    int q;

public:
    parameter(int p_, int q_) : p(p_), q(q_) {
        assert(p <= maxP);
        assert(p >= minP);
        assert(q >= 0);
        assert(p+q <= maxPplusQ);
    }

    int getP() const {return p;};
    int getQ() const {return q;};
};


string getFilePrefix(const parameter& par) {
    stringstream s;
    s << "../data/";
    s << par.getP();
    s << "_";
    s << par.getQ();
    s << "_";
    return s.str();
}

vector<parameter> readParameters() {
    vector<parameter> parameters;
    ifstream f("param.txt");
    string line;
    while (getline(f, line)) {
        istringstream iss(line);
        int precision;
        int maxRegisterValue;
        iss >> precision;
        iss >> maxRegisterValue;
        parameters.push_back(parameter(precision, maxRegisterValue));
    }
    return parameters;
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

int main(int argc, char* argv[])
{

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

    const vector<parameter> parameters = readParameters();

    for (size_t parameterIdx = 0; parameterIdx < parameters.size(); ++parameterIdx) {

        const parameter& par = parameters[parameterIdx];

        const unsigned char p = par.getP();
        const unsigned char q = par.getQ();

        const size_t resultBufferSize = seeds.size()*cardinalities.size();
        vector<double> flajoletSmallRangeEstimates(resultBufferSize);
        vector<double> flajoletMidRangeEstimates(resultBufferSize);
        vector<double> flajoletMidRangeEstimatesCorrected(resultBufferSize);
        vector<int> numSmallCorrectionIterations(resultBufferSize);
        vector<int> numLargeCorrectionIterations(resultBufferSize);
        vector<double> flajoletEstimates(resultBufferSize);
        vector<double> maxLikelihoodEstimates(resultBufferSize);
        vector<int> outerLoopIterationsCount(resultBufferSize);
        vector<int> innerLoop1IterationsCount(resultBufferSize);
        vector<int> innerLoop2IterationsCount(resultBufferSize);
        vector<int> logEvaluationCount(resultBufferSize);

        vector<double> strongLowerBoundEstimates(resultBufferSize);
        vector<double> weakLowerBoundEstimates(resultBufferSize);
        vector<double> strongUpperBoundEstimates(resultBufferSize);
        vector<double> weakUpperBoundEstimates(resultBufferSize);

        vector<int> kMinSum(cardinalities.size());
        vector<int> kMaxSum(cardinalities.size());

        const MaxLikelihoodEstimator maxLikelihoodEstimator(p,q);
        const CorrectedRawEstimator correctedRawEstimator(p,q);

        #pragma omp parallel for
        for(size_t seedCounter = 0; seedCounter < seeds.size(); ++seedCounter) {

            cout << getFileName(seeds[seedCounter], p) << endl;

            ifstream dataFile(getFileName(seeds[seedCounter], p));

            for(size_t cardinalityIdx = 0; cardinalityIdx < cardinalities.size(); ++cardinalityIdx) {

                string line;
                getline(dataFile, line);
                vector<int> tmpC = readCounts(line);
                vector<int> c(q+2);
                for (int i = 0; i < tmpC.size(); ++i) {
                    c[min(i, q+1)] += tmpC[i];
                }

                size_t resultPos = cardinalityIdx*seeds.size()+seedCounter;

                int kMin;
                int kMax;

                flajoletSmallRangeEstimates[resultPos] = flajoletSmallRangeEstimate(c);
                flajoletMidRangeEstimates[resultPos] = flajoletRawEstimate(c);
                flajoletMidRangeEstimatesCorrected[resultPos] = correctedRawEstimator.estimate(c, numSmallCorrectionIterations[resultPos], numLargeCorrectionIterations[resultPos]);
                flajoletEstimates[resultPos] = flajoletEstimate(c);
                maxLikelihoodEstimates[resultPos] = maxLikelihoodEstimator.estimate(c, outerLoopIterationsCount[resultPos], innerLoop1IterationsCount[resultPos], innerLoop2IterationsCount[resultPos], logEvaluationCount[resultPos], kMin, kMax);
                weakLowerBoundEstimates[resultPos] = weakLowerBoundEstimate(c);
                strongLowerBoundEstimates[resultPos] = strongLowerBoundEstimate(c);
                weakUpperBoundEstimates[resultPos] = weakUpperBoundEstimate(c);
                strongUpperBoundEstimates[resultPos] = strongUpperBoundEstimate(c);

                #pragma omp atomic
                kMinSum[cardinalityIdx] += kMin;

                #pragma omp atomic
                kMaxSum[cardinalityIdx] += kMax;
            }
        }

        vector<double> kMinAvg(cardinalities.size());
        vector<double> kMaxAvg(cardinalities.size());

        for(size_t cardinalityIdx = 0; cardinalityIdx < cardinalities.size(); ++cardinalityIdx) {
            kMinAvg[cardinalityIdx] = kMinSum[cardinalityIdx] / (double) seeds.size();
            kMaxAvg[cardinalityIdx] = kMaxSum[cardinalityIdx] / (double) seeds.size();
        }

        const string filePrefix = getFilePrefix(par);

        printToFile(filePrefix + "max_likelihood_estimates.dat", &maxLikelihoodEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "flajolet_small_range_estimates.dat", &flajoletSmallRangeEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "flajolet_mid_range_estimates.dat", &flajoletMidRangeEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "flajolet_mid_range_estimates_corrected.dat", &flajoletMidRangeEstimatesCorrected[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "num_small_correction_iterations.dat", &numSmallCorrectionIterations[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "num_large_correction_iterations.dat", &numLargeCorrectionIterations[0], seeds.size(), cardinalities.size());

        printToFile(filePrefix + "flajolet_estimates.dat", &flajoletEstimates[0], seeds.size(), cardinalities.size());

        printToFile(filePrefix + "outer_loop_iterations_count.dat", &outerLoopIterationsCount[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "inner_loop_1_iterations_count.dat", &innerLoop1IterationsCount[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "inner_loop_2_iterations_count.dat", &innerLoop2IterationsCount[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "log_evaluation_count.dat", &logEvaluationCount[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "kMin.dat", &kMinAvg[0], 1, cardinalities.size());
        printToFile(filePrefix + "kMax.dat", &kMaxAvg[0], 1, cardinalities.size());

        printToFile(filePrefix + "weak_lower_bound_estimates.dat", &weakLowerBoundEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "strong_lower_bound_estimates.dat", &strongLowerBoundEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "weak_upper_bound_estimates.dat", &weakUpperBoundEstimates[0], seeds.size(), cardinalities.size());
        printToFile(filePrefix + "strong_upper_bound_estimates.dat", &strongUpperBoundEstimates[0], seeds.size(), cardinalities.size());
    }
}
