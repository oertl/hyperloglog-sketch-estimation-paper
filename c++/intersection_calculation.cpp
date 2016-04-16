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
const string jointCardFileName = intersectDataPathName + "joint_cardinalities.dat";
const string resultsFileName = intersectDataPathName + "results.dat";


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
	
	resultsFile << "trueCardA;";
	resultsFile << "trueCardB;";
	resultsFile << "trueCardX;";
	
	resultsFile << "jaccardIndex;";
	resultsFile << "logRatio;";
	
	resultsFile << "inclExclMinA;";
	resultsFile << "inclExclM3sA;";
	resultsFile << "inclExclM2sA;";
	resultsFile << "inclExclM1sA;";
	resultsFile << "inclExclMedA;";
	resultsFile << "inclExclP1sA;";
	resultsFile << "inclExclP2sA;";
	resultsFile << "inclExclP3sA;";
	resultsFile << "inclExclMaxA;";
	
	resultsFile << "inclExclMinB;";
	resultsFile << "inclExclM3sB;";
	resultsFile << "inclExclM2sB;";
	resultsFile << "inclExclM1sB;";
	resultsFile << "inclExclMedB;";
	resultsFile << "inclExclP1sB;";
	resultsFile << "inclExclP2sB;";
	resultsFile << "inclExclP3sB;";
	resultsFile << "inclExclMaxB;";

	resultsFile << "inclExclMinX;";
	resultsFile << "inclExclM3sX;";
	resultsFile << "inclExclM2sX;";
	resultsFile << "inclExclM1sX;";
	resultsFile << "inclExclMedX;";
	resultsFile << "inclExclP1sX;";
	resultsFile << "inclExclP2sX;";
	resultsFile << "inclExclP3sX;";
	resultsFile << "inclExclMaxX;";
	
	resultsFile << "maxLikeMinA;";
	resultsFile << "maxLikeM3sA;";
	resultsFile << "maxLikeM2sA;";
	resultsFile << "maxLikeM1sA;";
	resultsFile << "maxLikeMedA;";
	resultsFile << "maxLikeP1sA;";
	resultsFile << "maxLikeP2sA;";
	resultsFile << "maxLikeP3sA;";
	resultsFile << "maxLikeMaxA;";
	
	resultsFile << "maxLikeMinB;";
	resultsFile << "maxLikeM3sB;";
	resultsFile << "maxLikeM2sB;";
	resultsFile << "maxLikeM1sB;";
	resultsFile << "maxLikeMedB;";
	resultsFile << "maxLikeP1sB;";
	resultsFile << "maxLikeP2sB;";
	resultsFile << "maxLikeP3sB;";
	resultsFile << "maxLikeMaxB;";

	resultsFile << "maxLikeMinX;";
	resultsFile << "maxLikeM3sX;";
	resultsFile << "maxLikeM2sX;";
	resultsFile << "maxLikeM1sX;";
	resultsFile << "maxLikeMedX;";
	resultsFile << "maxLikeP1sX;";
	resultsFile << "maxLikeP2sX;";
	resultsFile << "maxLikeP3sX;";
	resultsFile << "maxLikeMaxX;";
	
	resultsFile << "maxNumIterationsReachedCount;";
	resultsFile << "iterationAbortedCount" << endl;
	
	
	
	
	for(string fileName : jointCardFileNames) {
		
		const long trueCardA = stol(fileName.substr(0, 11));
		const long trueCardB = stol(fileName.substr(12, 23));
		const long trueCardX = stol(fileName.substr(24, 35));
		int p = stoi(fileName.substr(36, 38));
		int q = stoi(fileName.substr(39, 41));
		
		const double jaccardIndex = trueCardX/static_cast<double>(trueCardA+trueCardB+trueCardX);
		const double logRatio = std::log10(trueCardA) - std::log10(trueCardB);
		
		if (jaccardIndex < 1e-3 || jaccardIndex > 0.1 || std::fabs(logRatio) > 1) continue;

		cout << trueCardA << " " << trueCardB << " " << trueCardX << " " << p << " " << q << endl;
		
		ifstream statisticFile(intersectDataPathName + fileName);
		string line;
		
		std::vector<double> inExclEstimatedCardA;
		std::vector<double> inExclEstimatedCardB;
		std::vector<double> inExclEstimatedCardX;
		std::vector<double> maxLikeEstimatedCardA;
		std::vector<double> maxLikeEstimatedCardB;
		std::vector<double> maxLikeEstimatedCardX;
		int maxNumIterationsReachedCount = 0;
		int iterationAbortedCount = 0;
		
		while (getline(statisticFile, line))
		{
			const vector<int> jointStatistic = readCounts(line);
			assert(jointStatistic.size() == 5*(q+2));
			
			{
				double estCardA = 0.;
				double estCardB = 0.;
				double estCardX = 0.;
				inclusionExclusionTwoHyperLogLogEstimation(jointStatistic, estCardA, estCardB, estCardX);
				inExclEstimatedCardA.push_back(estCardA/trueCardA-1.);
				inExclEstimatedCardB.push_back(estCardB/trueCardB-1.);
				inExclEstimatedCardX.push_back(estCardX/trueCardX-1.);
			}
			
			{
				double estCardA = 0.;
				double estCardB = 0.;
				double estCardX = 0.;
				bool maxNumIterationsReached;
				bool iterationAborted;
				maxLikelihoodTwoHyperLogLogEstimation2(jointStatistic, estCardA, estCardB, estCardX, maxNumIterationsReached, iterationAborted);
				maxLikeEstimatedCardA.push_back(estCardA/trueCardA-1.);
				maxLikeEstimatedCardB.push_back(estCardB/trueCardB-1.);
				maxLikeEstimatedCardX.push_back(estCardX/trueCardX-1.);
				if (maxNumIterationsReached) maxNumIterationsReachedCount+=1;
				if (iterationAborted) iterationAbortedCount += 1;
			}
			
		}
		
		std::sort(inExclEstimatedCardA.begin(), inExclEstimatedCardA.end());
		std::sort(inExclEstimatedCardB.begin(), inExclEstimatedCardB.end());
		std::sort(inExclEstimatedCardX.begin(), inExclEstimatedCardX.end());
		std::sort(maxLikeEstimatedCardA.begin(), maxLikeEstimatedCardA.end());
		std::sort(maxLikeEstimatedCardB.begin(), maxLikeEstimatedCardB.end());
		std::sort(maxLikeEstimatedCardX.begin(), maxLikeEstimatedCardX.end());
		
		resultsFile << trueCardA << ";";
		resultsFile << trueCardB << ";";
		resultsFile << trueCardX << ";";
		
		resultsFile << jaccardIndex << ";";
		resultsFile << logRatio << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardA[0], 1, inExclEstimatedCardA.size(), 1.0) << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardB[0], 1, inExclEstimatedCardB.size(), 1.0) << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&inExclEstimatedCardX[0], 1, inExclEstimatedCardX.size(), 1.0) << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardA[0], 1, maxLikeEstimatedCardA.size(), 1.0) << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardB[0], 1, maxLikeEstimatedCardB.size(), 1.0) << ";";
		
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), 0.0) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pm3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pm2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pm1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), 0.5) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pp1s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pp2s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), pp3s) << ";";
		resultsFile << gsl_stats_quantile_from_sorted_data(&maxLikeEstimatedCardX[0], 1, maxLikeEstimatedCardX.size(), 1.0) << ";";
		
		resultsFile << maxNumIterationsReachedCount << ";";
		resultsFile << iterationAbortedCount << endl;
	}
}
