//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#include "cardinality_estimation.hpp"
#include "two_hyperloglog_statistic.hpp"
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
#include <gsl/gsl_statistics.h>

using namespace std;

const int minP = 4;
const int maxP = 22;
const int maxPplusQ = 64;

const double pm3s = (1.-erf(3./sqrt(2.)))/2.;
const double pp3s = (1.+erf(3./sqrt(2.)))/2.;
const double pm2s = (1.-erf(2./sqrt(2.)))/2.;
const double pp2s = (1.+erf(2./sqrt(2.)))/2.;
const double pm1s = (1.-erf(1./sqrt(2.)))/2.;
const double pp1s = (1.+erf(1./sqrt(2.)))/2.;
const double median = 0.5;


const string intersectDataPathName = "../data/hll_joint/";
const string paperPathName = "../paper/";
const string jointCardFileName = intersectDataPathName + "joint_cardinalities.dat";
const string resultsFileName = paperPathName + "intersection.csv";


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

void printHeaders(ofstream& os, string prefix, string postfix) {
    os << prefix << "Min" << postfix << ",";
    os << prefix << "M3s" << postfix << ",";
    os << prefix << "M2s" << postfix << ",";
    os << prefix << "M1s" << postfix << ",";
    os << prefix << "Med" << postfix << ",";
    os << prefix << "P1s" << postfix << ",";
    os << prefix << "P2s" << postfix << ",";
    os << prefix << "P3s" << postfix << ",";
    os << prefix << "Max" << postfix << ",";
    os << prefix << "Mean" << postfix << ",";
    os << prefix << "StdDev" << postfix << ",";
}

void calculateAndPrintStatistics(ofstream& os, std::vector<double> data) {
    std::sort(data.begin(), data.end());
    os << gsl_stats_min(&data[0], 1, data.size()) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pm3s) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pm2s) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pm1s) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), 0.5) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pp1s) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pp2s) << ",";
    os << gsl_stats_quantile_from_sorted_data(&data[0], 1, data.size(), pp3s) << ",";
    os << gsl_stats_max(&data[0], 1, data.size()) << ",";
    const double mean = gsl_stats_mean(&data[0], 1, data.size());
    os << mean << ",";
    os << gsl_stats_sd_m(&data[0], 1, data.size(), mean) << ",";
}

int main(int argc, char* argv[])
{

    vector<string> jointCardFileNames;
    {
        ifstream jointCardFileNamesFile(jointCardFileName);

        string line;
        while (getline(jointCardFileNamesFile, line))
        {
            jointCardFileNames.push_back(line);
        }
    }

    cout << "jointCardFileNames size "  << jointCardFileNames.size() << endl;

    ofstream resultsFile(resultsFileName);

    resultsFile << "trueCardA,";
    resultsFile << "trueCardB,";
    resultsFile << "trueCardX,";

    resultsFile << "trueJaccardIdx,";
    resultsFile << "trueLogRatio,";

    printHeaders(resultsFile, "inclExcl", "A");
    printHeaders(resultsFile, "inclExcl", "B");
    printHeaders(resultsFile, "inclExcl", "X");
    printHeaders(resultsFile, "inclExcl", "JaccardIdx");
    printHeaders(resultsFile, "maxLike", "A");
    printHeaders(resultsFile, "maxLike", "B");
    printHeaders(resultsFile, "maxLike", "X");
    printHeaders(resultsFile, "maxLike", "JaccardIdx");

    resultsFile << "maxNumIterationsReachedCount,";
    resultsFile << "avgNumIterations,";
    resultsFile << "maxNumIterations,";
    resultsFile << "iterationAbortedCount" << endl;

    for(const string& fileName : jointCardFileNames) {

        const long trueCardA = stol(fileName.substr(0, 11));
        const long trueCardB = stol(fileName.substr(12, 23));
        const long trueCardX = stol(fileName.substr(24, 35));
        int p = stoi(fileName.substr(36, 38));
        int q = stoi(fileName.substr(39, 41));

        const double jaccardIndex = trueCardX/static_cast<double>(trueCardA+trueCardB+trueCardX);
        const double logRatio = std::log10(trueCardA) - std::log10(trueCardB);

        if (jaccardIndex < 1e-3 || jaccardIndex > 0.1 || std::fabs(logRatio) > 1) continue; // TODO filter interesting cardinalities

        cout << trueCardA << " " << trueCardB << " " << trueCardX << " " << p << " " << q << endl;

        ifstream statisticFile(intersectDataPathName + fileName);
        string line;

        std::vector<double> inExclEstimatedCardA;
        std::vector<double> inExclEstimatedCardB;
        std::vector<double> inExclEstimatedCardX;
        std::vector<double> inExclJaccardIdx;
        std::vector<double> maxLikeEstimatedCardA;
        std::vector<double> maxLikeEstimatedCardB;
        std::vector<double> maxLikeEstimatedCardX;
        std::vector<double> maxLikeJaccardIdx;
        int maxNumIterationsReachedCount = 0;
        int maxNumIterations = 0;
        int iterationAbortedCount = 0;
        long numIterationsTotal = 0;
        int size = 0;

        while (getline(statisticFile, line))
        {
            TwoHyperLogLogStatistic jointStatistic = TwoHyperLogLogStatistic::fromString(line);
            assert(jointStatistic.getQ() == q);

            {
                double estCardA = 0.;
                double estCardB = 0.;
                double estCardX = 0.;
                inclusionExclusionTwoHyperLogLogEstimation(jointStatistic, estCardA, estCardB, estCardX);
                inExclEstimatedCardA.push_back(estCardA/trueCardA-1.);
                inExclEstimatedCardB.push_back(estCardB/trueCardB-1.);
                inExclEstimatedCardX.push_back(estCardX/trueCardX-1.);
                inExclJaccardIdx.push_back((estCardX/(estCardA+estCardB+estCardX))/jaccardIndex-1.);
            }

            {
                double estCardA = 0.;
                double estCardB = 0.;
                double estCardX = 0.;
                bool maxNumIterationsReached = false;
                bool iterationAborted = false;
                int numIterations = 0;
                //analyticalJointHyperLogLogEstimator(jointStatistic, estCardA, estCardB, estCardX);
                maxLikelihoodTwoHyperLogLogEstimation(jointStatistic, estCardA, estCardB, estCardX, maxNumIterationsReached, iterationAborted, numIterations);
                maxLikeEstimatedCardA.push_back(estCardA/trueCardA-1.);
                maxLikeEstimatedCardB.push_back(estCardB/trueCardB-1.);
                maxLikeEstimatedCardX.push_back(estCardX/trueCardX-1.);
                maxLikeJaccardIdx.push_back((estCardX/(estCardA+estCardB+estCardX))/jaccardIndex-1.);
                if (maxNumIterationsReached) maxNumIterationsReachedCount+=1;
                if (iterationAborted) iterationAbortedCount += 1;
                numIterationsTotal += numIterations;
                maxNumIterations = std::max(maxNumIterations, numIterations);
            }

            size += 1;

        }

        resultsFile << trueCardA << ",";
        resultsFile << trueCardB << ",";
        resultsFile << trueCardX << ",";

        resultsFile << jaccardIndex << ",";
        resultsFile << logRatio << ",";

        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardA);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardB);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardX);
        calculateAndPrintStatistics(resultsFile, inExclJaccardIdx);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardA);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardB);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardX);
        calculateAndPrintStatistics(resultsFile, maxLikeJaccardIdx);

        resultsFile << maxNumIterationsReachedCount << ",";
        resultsFile << (numIterationsTotal/(double)size) << ",";
        resultsFile << maxNumIterations << ",";
        resultsFile << iterationAbortedCount << endl;
    }
}
