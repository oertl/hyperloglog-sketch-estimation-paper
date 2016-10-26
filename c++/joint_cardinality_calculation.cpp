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
    os << prefix << "RMSE" << postfix << ",";
}

void calculateAndPrintStatistics(ofstream& os, std::vector<double> data, double& stdev, double& rmse) {
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
    stdev = gsl_stats_sd_with_fixed_mean(&data[0], 1, data.size(), mean);
    rmse = gsl_stats_sd_with_fixed_mean(&data[0], 1, data.size(), 0);
    os << mean << ",";
    os << stdev << ",";
    os << rmse << ",";
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

    resultsFile << "avgNumEvaluations,";
    resultsFile << "avgNumIterations,";
    resultsFile << "maxNumIterations" << endl;

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
        int maxNumIterations = 0;
        long numIterationsTotal = 0;
        long numEvaluationsTotal = 0;
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
                inExclEstimatedCardAX.push_back((estCardA+estCardX)/(trueCardA+trueCardX)-1.);
                inExclEstimatedCardBX.push_back((estCardB+estCardX)/(trueCardB+trueCardX)-1.);
                inExclEstimatedCardABX.push_back((estCardA+estCardB+estCardX)/(trueCardA+trueCardB+trueCardX)-1.);
                inExclJaccardIdx.push_back((estCardX/(estCardA+estCardB+estCardX))/jaccardIndex-1.);
            }

            {
                double estCardA = 0.;
                double estCardB = 0.;
                double estCardX = 0.;
                int numIterations = 0;
                int numEvaluations = 0;
                maxLikelihoodTwoHyperLogLogEstimation(jointStatistic, estCardA, estCardB, estCardX, numIterations, numEvaluations);
                maxLikeEstimatedCardA.push_back(estCardA/trueCardA-1.);
                maxLikeEstimatedCardB.push_back(estCardB/trueCardB-1.);
                maxLikeEstimatedCardX.push_back(estCardX/trueCardX-1.);
                maxLikeJaccardIdx.push_back((estCardX/(estCardA+estCardB+estCardX))/jaccardIndex-1.);
                maxLikeEstimatedCardAX.push_back((estCardA+estCardX)/(trueCardA+trueCardX)-1.);
                maxLikeEstimatedCardBX.push_back((estCardB+estCardX)/(trueCardB+trueCardX)-1.);
                maxLikeEstimatedCardABX.push_back((estCardA+estCardB+estCardX)/(trueCardA+trueCardB+trueCardX)-1.);
                numIterationsTotal += numIterations;
                maxNumIterations = std::max(maxNumIterations, numIterations);
                numEvaluationsTotal += numEvaluations;
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

        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardA, inExclStdevA, inExclRmseA);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardB, inExclStdevB, inExclRmseB);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardX, inExclStdevX, inExclRmseX);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardAX, inExclStdevAX, inExclRmseAX);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardBX, inExclStdevBX, inExclRmseBX);
        calculateAndPrintStatistics(resultsFile, inExclEstimatedCardABX, inExclStdevABX, inExclRmseABX);
        calculateAndPrintStatistics(resultsFile, inExclJaccardIdx, inExclStdevJaccard, inExclRmseJaccard);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardA, maxLikeStdevA, maxLikeRmseA);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardB, maxLikeStdevB, maxLikeRmseB);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardX, maxLikeStdevX, maxLikeRmseX);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardAX, maxLikeStdevAX, maxLikeRmseAX);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardBX, maxLikeStdevBX, maxLikeRmseBX);
        calculateAndPrintStatistics(resultsFile, maxLikeEstimatedCardABX, maxLikeStdevABX, maxLikeRmseABX);
        calculateAndPrintStatistics(resultsFile, maxLikeJaccardIdx, maxLikeStdevJaccard, maxLikeRmseJaccard);

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

        resultsFile << (numEvaluationsTotal/(double)size) << ",";
        resultsFile << (numIterationsTotal/(double)size) << ",";
        resultsFile << maxNumIterations << endl;

    }
    cout << "number of evaluations = " << evaluationCounter << endl;
}
