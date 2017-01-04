//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#define CARDINALITY_ESTIMATION_USE_DLIB

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

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

using namespace std;
using namespace boost::accumulators;

const int minP = 4;
const int maxP = 22;
const int maxPplusQ = 64;

const int analyzeP = 20;

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
const string resultsFileName = paperPathName + "joint_cardinality_calculation.csv";


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
    os << prefix << "Max" << postfix << ",";
    os << prefix << "Mean" << postfix << ",";
    os << prefix << "StdDev" << postfix << ",";
    os << prefix << "RMSE" << postfix << ",";
}

void calculateAndPrintStatistics(ofstream& os, double trueValue, std::vector<double> data, double& relativeStdev, double& relativeRmse) {

    accumulator_set<double, stats<tag::min, tag::max, tag::sum_kahan> > acc;

    for (auto d : data) {
        acc(d);
    }

    const double mean = extract_result<tag::sum_kahan>(acc)/data.size();
    const double min = extract_result<tag::min>(acc);
    const double max = extract_result<tag::max>(acc);

    accumulator_set<double, stats<tag::sum_kahan> > acc_rmse;
    accumulator_set<double, stats<tag::sum_kahan> > acc_stdev;

    for (auto d : data) {
        acc_rmse(std::pow(d - trueValue, 2));
        acc_stdev(std::pow(d - mean, 2));
    }

    const double rmse = std::sqrt(extract_result<tag::sum_kahan>(acc_rmse)/data.size());
    const double stdev = std::sqrt(extract_result<tag::sum_kahan>(acc_stdev)/data.size());

    const double relativeMean = mean/trueValue - 1.;
    const double relativeMin = min/trueValue - 1.;
    const double relativeMax = max/trueValue - 1.;
    relativeRmse = rmse/trueValue;
    relativeStdev = stdev/trueValue;

    os << relativeMin << ",";
    os << relativeMax << ",";
    os << relativeMean << ",";
    os << relativeStdev << ",";
    os << relativeRmse << ",";
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
    resultsFile << "trueRatio,";

    printHeaders(resultsFile, "inclExcl", "A");
    printHeaders(resultsFile, "inclExcl", "B");
    printHeaders(resultsFile, "inclExcl", "X");
    printHeaders(resultsFile, "inclExcl", "AX");
    printHeaders(resultsFile, "inclExcl", "BX");
    printHeaders(resultsFile, "inclExcl", "ABX");
    printHeaders(resultsFile, "inclExcl", "JaccardIdx");
    printHeaders(resultsFile, "maxLike", "A");
    printHeaders(resultsFile, "maxLike", "B");
    printHeaders(resultsFile, "maxLike", "X");
    printHeaders(resultsFile, "maxLike", "AX");
    printHeaders(resultsFile, "maxLike", "BX");
    printHeaders(resultsFile, "maxLike", "ABX");
    printHeaders(resultsFile, "maxLike", "JaccardIdx");

    resultsFile << "improvementStdevA,";
    resultsFile << "improvementRmseA,";
    resultsFile << "improvementStdevB,";
    resultsFile << "improvementRmseB,";
    resultsFile << "improvementStdevX,";
    resultsFile << "improvementRmseX,";
    resultsFile << "improvementStdevAX,";
    resultsFile << "improvementRmseAX,";
    resultsFile << "improvementStdevBX,";
    resultsFile << "improvementRmseBX,";
    resultsFile << "improvementStdevABX,";
    resultsFile << "improvementRmseABX,";
    resultsFile << "improvementStdevJaccard,";
    resultsFile << "improvementRmseJaccard,";

    resultsFile << "avgNumFunctionEvaluations,";
    resultsFile << "avgNumGradientEvaluations" << endl;

    int evaluationCounter = 0;

    for(const string& fileName : jointCardFileNames) {

        long trueCardA = stol(fileName.substr(0, 11));
        long trueCardB = stol(fileName.substr(12, 23));

        const long trueCardX = stol(fileName.substr(24, 35));
        int p = stoi(fileName.substr(36, 38));
        int q = stoi(fileName.substr(39, 41));

        const double jaccardIndex = trueCardX/static_cast<double>(trueCardA+trueCardB+trueCardX);
        double logRatio = std::log10(trueCardA) - std::log10(trueCardB);

        if (analyzeP != p || jaccardIndex < 1e-3 || jaccardIndex > 0.2 || std::fabs(logRatio) > 2) continue; // TODO filter cardinality combinations for table

        cout << trueCardA << " " << trueCardB << " " << trueCardX << " " << p << " " << q << endl;

        ifstream statisticFile(intersectDataPathName + fileName);
        string line;

        std::vector<double> inExclEstimatedCardA;
        std::vector<double> inExclEstimatedCardB;
        std::vector<double> inExclEstimatedCardX;
        std::vector<double> inExclEstimatedCardAX;
        std::vector<double> inExclEstimatedCardBX;
        std::vector<double> inExclEstimatedCardABX;
        std::vector<double> inExclJaccardIdx;
        std::vector<double> maxLikeEstimatedCardA;
        std::vector<double> maxLikeEstimatedCardB;
        std::vector<double> maxLikeEstimatedCardX;
        std::vector<double> maxLikeEstimatedCardAX;
        std::vector<double> maxLikeEstimatedCardBX;
        std::vector<double> maxLikeEstimatedCardABX;
        std::vector<double> maxLikeJaccardIdx;
        size_t numFunctionEvaluationsTotal = 0;
        size_t numGradientEvaluationsTotal = 0;
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
                inExclEstimatedCardA.push_back(estCardA);
                inExclEstimatedCardB.push_back(estCardB);
                inExclEstimatedCardX.push_back(estCardX);
                inExclEstimatedCardAX.push_back(estCardA+estCardX);
                inExclEstimatedCardBX.push_back(estCardB+estCardX);
                inExclEstimatedCardABX.push_back(estCardA+estCardB+estCardX);
                inExclJaccardIdx.push_back(estCardX/(estCardA+estCardB+estCardX));
            }

            {
                double estCardA = 0.;
                double estCardB = 0.;
                double estCardX = 0.;
                size_t numFunctionEvaluations = 0;
                size_t numGradientEvaluations = 0;
                maxLikelihoodTwoHyperLogLogEstimation(jointStatistic, estCardA, estCardB, estCardX, numFunctionEvaluations, numGradientEvaluations);
                assert(numFunctionEvaluations == numGradientEvaluations); // TODO
                maxLikeEstimatedCardA.push_back(estCardA);
                maxLikeEstimatedCardB.push_back(estCardB);
                maxLikeEstimatedCardX.push_back(estCardX);
                maxLikeEstimatedCardAX.push_back(estCardA+estCardX);
                maxLikeEstimatedCardBX.push_back(estCardB+estCardX);
                maxLikeEstimatedCardABX.push_back(estCardA+estCardB+estCardX);
                maxLikeJaccardIdx.push_back(estCardX/(estCardA+estCardB+estCardX));
                numFunctionEvaluationsTotal += numFunctionEvaluations;
                numGradientEvaluationsTotal += numGradientEvaluations;

            }

            size += 1;
            evaluationCounter += 1;

        }

        if(logRatio < 0) {
            swap(inExclEstimatedCardA, inExclEstimatedCardB);
            swap(inExclEstimatedCardAX, inExclEstimatedCardBX);
            swap(maxLikeEstimatedCardA, maxLikeEstimatedCardB);
            swap(maxLikeEstimatedCardAX, maxLikeEstimatedCardBX);
            swap(trueCardA, trueCardB);
            logRatio = -logRatio;
        }
        const double ratio = static_cast<double>(trueCardA)/static_cast<double>(trueCardB);

        resultsFile << trueCardA << ",";
        resultsFile << trueCardB << ",";
        resultsFile << trueCardX << ",";

        resultsFile << jaccardIndex << ",";
        resultsFile << ratio << ",";

        double inExclRmseA, inExclStdevA, inExclRmseB, inExclStdevB, inExclRmseX, inExclStdevX, inExclRmseAX, inExclStdevAX, inExclRmseBX, inExclStdevBX, inExclRmseABX, inExclStdevABX;
        double maxLikeRmseA, maxLikeStdevA, maxLikeRmseB, maxLikeStdevB, maxLikeRmseX, maxLikeStdevX, maxLikeRmseAX, maxLikeStdevAX, maxLikeRmseBX, maxLikeStdevBX, maxLikeRmseABX, maxLikeStdevABX;

        double inExclRmseJaccard, inExclStdevJaccard;
        double maxLikeRmseJaccard, maxLikeStdevJaccard;

        calculateAndPrintStatistics(resultsFile, trueCardA, inExclEstimatedCardA, inExclStdevA, inExclRmseA);
        calculateAndPrintStatistics(resultsFile, trueCardB, inExclEstimatedCardB, inExclStdevB, inExclRmseB);
        calculateAndPrintStatistics(resultsFile, trueCardX, inExclEstimatedCardX, inExclStdevX, inExclRmseX);
        calculateAndPrintStatistics(resultsFile, trueCardA + trueCardX, inExclEstimatedCardAX, inExclStdevAX, inExclRmseAX);
        calculateAndPrintStatistics(resultsFile, trueCardB + trueCardX, inExclEstimatedCardBX, inExclStdevBX, inExclRmseBX);
        calculateAndPrintStatistics(resultsFile, trueCardA + trueCardB + trueCardX, inExclEstimatedCardABX, inExclStdevABX, inExclRmseABX);
        calculateAndPrintStatistics(resultsFile, jaccardIndex, inExclJaccardIdx, inExclStdevJaccard, inExclRmseJaccard);
        calculateAndPrintStatistics(resultsFile, trueCardA, maxLikeEstimatedCardA, maxLikeStdevA, maxLikeRmseA);
        calculateAndPrintStatistics(resultsFile, trueCardB, maxLikeEstimatedCardB, maxLikeStdevB, maxLikeRmseB);
        calculateAndPrintStatistics(resultsFile, trueCardX, maxLikeEstimatedCardX, maxLikeStdevX, maxLikeRmseX);
        calculateAndPrintStatistics(resultsFile, trueCardA + trueCardX, maxLikeEstimatedCardAX, maxLikeStdevAX, maxLikeRmseAX);
        calculateAndPrintStatistics(resultsFile, trueCardB + trueCardX, maxLikeEstimatedCardBX, maxLikeStdevBX, maxLikeRmseBX);
        calculateAndPrintStatistics(resultsFile, trueCardA + trueCardB + trueCardX, maxLikeEstimatedCardABX, maxLikeStdevABX, maxLikeRmseABX);
        calculateAndPrintStatistics(resultsFile, jaccardIndex, maxLikeJaccardIdx, maxLikeStdevJaccard, maxLikeRmseJaccard);

        resultsFile << inExclStdevA/maxLikeStdevA  << ",";
        resultsFile << inExclRmseA/maxLikeRmseA  << ",";
        resultsFile << inExclStdevB/maxLikeStdevB  << ",";
        resultsFile << inExclRmseB/maxLikeRmseB  << ",";
        resultsFile << inExclStdevX/maxLikeStdevX  << ",";
        resultsFile << inExclRmseX/maxLikeRmseX  << ",";
        resultsFile << inExclStdevAX/maxLikeStdevAX  << ",";
        resultsFile << inExclRmseAX/maxLikeRmseAX  << ",";
        resultsFile << inExclStdevBX/maxLikeStdevBX  << ",";
        resultsFile << inExclRmseBX/maxLikeRmseBX  << ",";
        resultsFile << inExclStdevABX/maxLikeStdevABX  << ",";
        resultsFile << inExclRmseABX/maxLikeRmseABX  << ",";
        resultsFile << inExclStdevJaccard/maxLikeStdevJaccard << ",";
        resultsFile << inExclRmseJaccard/maxLikeRmseJaccard << ",";

        resultsFile << (numFunctionEvaluationsTotal/(double)size) << ",";
        resultsFile << (numGradientEvaluationsTotal/(double)size) << endl;

    }
    cout << "number of evaluations = " << evaluationCounter << endl;
}
