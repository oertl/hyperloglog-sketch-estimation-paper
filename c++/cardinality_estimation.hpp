//#############################
//# Copyright 2016 Otmar Ertl #
//#############################

#ifndef _CARDINALITY_ESTIMATION_HPP_
#define _CARDINALITY_ESTIMATION_HPP_

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <functional>

#ifdef CARDINALITY_ESTIMATION_USE_CERES
#include "ceres/ceres.h"
#endif

#ifdef CARDINALITY_ESTIMATION_USE_DLIB
#include "dlib/optimization.h"
#endif

#include "two_hyperloglog_statistic.hpp"

int getPFromNumberRegisters(int m) {
    int p;
    std::frexp(m, &p);
    p -= 1;
    assert(m == 1 << p);
    return p;
}

int getPFromCounts(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    return getPFromNumberRegisters(m);
}

class MaxLikelihoodEstimator {
    const double eps = 1e-2;
    const int p;
    const int q;
    const int m;
    const double relativeErrorLimit;

public:

    MaxLikelihoodEstimator(const int p_, const int q_) :p(p_), q(q_), m(1 << p_), relativeErrorLimit(eps/(sqrt(m))) {}

    double estimate(const std::vector<int>& c, int& outerLoopIterationsCount, int& innerLoop1IterationsCount, int& innerLoop2IterationsCount, int& logEvaluationCount) const {

        outerLoopIterationsCount = 0;
        innerLoop1IterationsCount = 0;
        innerLoop2IterationsCount = 0;
        logEvaluationCount = 0;

        if (c[q+1] == m) return std::numeric_limits<double>::infinity();

        int kMin, kMax;

        for(kMin=0; c[kMin]==0; ++kMin);
        int kMinPrime = std::max(1, kMin);

        for(kMax=q+1; c[kMax]==0; --kMax);
        int kMaxPrime = std::min(q, kMax);

        double z = 0.;
        for (int k = kMaxPrime; k >= kMinPrime; --k) {
            z = 0.5*z + c[k];
        }
        z = ldexp(z, -kMinPrime);

        int cPrime = c[q+1];
        if (q >= 1) {
            cPrime += c[kMaxPrime];
        }

        double gPrevious = 0;

        double x;
        double a = z + c[0];
        int mPrime = m - c[0];
        {
            double b = z + ldexp(c[q+1], -q);
            if (b <= 1.5*a) {
                x = mPrime/(0.5*b+a);
            }
            else {
                x = (mPrime/b)*log1p(b/a);
                logEvaluationCount += 1;
            }
        }

        double deltaX = x;
        while(deltaX > x*relativeErrorLimit) {
            int kappaMinus1;
            frexp(x, &kappaMinus1);
            double xPrime = ldexp(x, -std::max(kMaxPrime+1, kappaMinus1+2));
            double xPrime2 = xPrime*xPrime;
            double h = xPrime - xPrime2/3 + (xPrime2*xPrime2)*(1./45. - xPrime2/472.5);
            for (int k = kappaMinus1; k >= kMaxPrime; --k) {
                double hPrime = 1. - h;
                h = (xPrime + h*hPrime)/(xPrime+hPrime);
                xPrime += xPrime;
                innerLoop1IterationsCount += 1;
            }
            double g = cPrime*h;
            for (int k = kMaxPrime-1; k >= kMinPrime; --k) {
                double hPrime = 1. - h;
                h = (xPrime + h*hPrime)/(xPrime+hPrime);
                xPrime += xPrime;
                g += c[k] * h;
                innerLoop2IterationsCount += 1;
            }
            g += x*a;
            if (gPrevious < g && g <= mPrime) {
                deltaX *= (g-mPrime)/(gPrevious-g);
            }
            else {
                deltaX = 0;
            }
            x += deltaX;
            gPrevious = g;
            outerLoopIterationsCount += 1;
        }
        return x*m;
    }

    double operator()(const std::vector<int>& c) const {

        int outerLoopIterationsCount;
        int innerLoop1IterationsCount;
        int innerLoop2IterationsCount;
        int logEvaluationCount;

        return estimate(c, outerLoopIterationsCount, innerLoop1IterationsCount, innerLoop2IterationsCount, logEvaluationCount);
    }
};

double flajoletSmallRangeEstimate(const std::vector<int>& c, int m) {
    return m*std::log(m/(double)(c[0]));
}

double flajoletSmallRangeEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    return flajoletSmallRangeEstimate(c, m);
}

double getAlpha(int m) {
    //return 1./(2.*std::log(2));
    assert(m >= 16);
    double alpha;
    if (m == 16) {
        alpha = 0.673;
    }
    else if (m == 32) {
        alpha = 0.697;
    }
    else if (m == 64) {
        alpha = 0.709;
    }
    else {
        alpha = 0.7213/(1. + 1.079/m);
    }
    return alpha;
}

double flajoletRawEstimate(const std::vector<int>& c, int m) {
    double z = 0;
    for (int i = c.size()-1; i >= 0; --i) {
        z += ldexp(c[i], -i);
    }
    double alpha = getAlpha(m);
    return (alpha*m)*m/z;
}

class ImprovedRawEstimator {
    const int p;
    const int q;
    const int m = 1 << p;
    const double m2alpha = (m/(2.*std::log(2)))*m;
    const std::vector<double> tauValues = initTau(m);
    const std::vector<double> sigmaValues = initSigma(m);

    static double sigma(double x, int& numIterations) {
        numIterations = 0;
        if (x == 1.) return std::numeric_limits<double>::infinity();
        double zPrime;
        double y = 1;
        double z = x;
        do {
            numIterations += 1;
            x *= x;
            zPrime = z;
            z += x * y;
            y += y;
        } while(zPrime != z);
        return z;
    }

    static double tau(double x, int& numIterations) {
        numIterations = 0;
        if (x == 0. || x == 1.) return 0.;
        double zPrime;
        double y = 1.0;
        double z = 0;
        do {
            numIterations += 1;
            x = std::sqrt(x);
            zPrime = z;
            y *= 0.5;
            z += (1 - x)*x*y;
        } while(zPrime != z);
        return z;
    }

    static std::vector<double> initSigma(int m) {
        int numIterations;
        std::vector<double> result(m+1);
        for (int c = 0; c <= m; ++c) {
            result[c] = m * sigma(static_cast<double>(c)/static_cast<double>(m), numIterations);
        }
        return result;
    }

    static std::vector<double> initTau(int m) {
        int numIterations;
        std::vector<double> result(m+1);
        for (int c = 0; c <= m; ++c) {
            result[c] = m * tau(static_cast<double>(m-c)/static_cast<double>(m), numIterations);
        }
        return result;
    }

public:

    ImprovedRawEstimator(const int p_, const int q_) : p(p_), q(q_) {}

    double estimate_on_demand(const std::vector<int>& c, int& numSmallCorrectionIterations, int& numLargeCorrectionIterations) const {

        numSmallCorrectionIterations = 0;
        numLargeCorrectionIterations = 0;

        double z = m * tau(static_cast<double>(m-c[q+1])/static_cast<double>(m), numLargeCorrectionIterations);
        for (int k = q; k >= 1; --k) {
            z += c[k];
            z *= 0.5;
        }
        z += m * sigma(static_cast<double>(c[0])/static_cast<double>(m), numSmallCorrectionIterations);
        return m2alpha/z;
    }

    double estimate_precalculated(const std::vector<int>& c) const {

        double z = tauValues[c[q+1]];
        for (int k = q; k >= 1; --k) {
            z += c[k];
            z *= 0.5;
        }
        z += sigmaValues[c[0]];
        return m2alpha/z;
    }

};

double flajoletRawEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    return flajoletRawEstimate(c, m);
}

double flajoletEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    int q = c.size()-2;
    int p = getPFromNumberRegisters(m);
    if (p + q != 32) {
        return 0.; // original estimate is only designed for 32 bit hash values
    }

    double e  = flajoletRawEstimate(c, m);
    double eStar;
    if (e <= 5./2.*m) {
        int v = c[0];
        if (v != 0) {
            eStar = m * std::log(static_cast<double>(m)/v);
        }
        else {
            eStar = e;
        }
    }
    else if (e <= ldexp(1./30., 32)) {
        eStar = e;
    }
    else {
        if (e <= ldexp(1., 32)) { // is not part of original algorithm
            eStar = -std::pow(2., 32)*std::log1p(-ldexp(e, -32));
        }
        else {
            eStar = std::numeric_limits<double>::infinity();
        }
    }
    return eStar;
}

double strongLowerBoundEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    int q = c.size()-2;
    if (m == c[0]) return 0.;
    if (m == c[q+1]) return std::numeric_limits<double>::infinity();

    double s = 0.;
    for (int i = q; i >= 1; --i) {
        s += std::ldexp(c[i], -i);
    }

    double a = s + c[0];
    double b = s + std::ldexp(c[q+1], -q);
    int numNonZeroRegisters = m-c[0];
    double x = (numNonZeroRegisters/b)*log1p(b/a);
    return m*x;
}

double weakLowerBoundEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    int q = c.size()-2;
    if (m == c[0]) return 0.;
    if (m == c[q+1]) return std::numeric_limits<double>::infinity();

    double s = 0.;
    for (int i = q; i >= 1; --i) {
        s += std::ldexp(c[i], -i);
    }

    double a = s + c[0];
    double b = s + ldexp(c[q+1], -q);
    int numNonZeroRegisters = m-c[0];
    double x = numNonZeroRegisters/(0.5*b+a);
    return m*x;
}

double weakUpperBoundEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    int q = c.size()-2;
    if (m == c[0]) return 0.;
    if (m == c[q+1]) return std::numeric_limits<double>::infinity();

    double z = 0;
    for (int i = q; i >= 0; --i) {
        z += std::ldexp(c[i], -i);
    }
    return ((m-c[0])/z)*m;
}

double strongUpperBoundEstimate(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    int q = c.size()-2;
    if (m == c[0]) return 0.;
    if (m == c[q+1]) return std::numeric_limits<double>::infinity();
    int kMax;
    for(kMax=q+1; c[kMax]==0; --kMax);

    int kPrime = std::min(kMax, q);

    double z = 0;
    for (int i = q; i >= 0; --i) {
        z += std::ldexp(c[i], -i);
    }
    return std::ldexp(m*std::log1p((m-c[0])/(ldexp(z, kPrime))), kPrime);
}

void newTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX) {

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    const double cardinalityAX = estimator(jointStatistic.get1Counts());
    const double cardinalityBX = estimator(jointStatistic.get2Counts());
    const double cardinalityABX = estimator(jointStatistic.getMaxCounts());


    std::vector<int> countsAXBhalf(jointStatistic.getQ() + 1);
    std::vector<int> countsBXAhalf(jointStatistic.getQ() + 1);
    int sumAXBhalf= jointStatistic.getNumRegisters();
    int sumBXAhalf= jointStatistic.getNumRegisters();
    for (int q = 0; q < jointStatistic.getQ(); ++q) {
        countsAXBhalf[q] = jointStatistic.getLarger1Count(q) + jointStatistic.getEqualCount(q) + jointStatistic.getLarger2Count(q+1);
        countsBXAhalf[q] = jointStatistic.getLarger2Count(q) + jointStatistic.getEqualCount(q) + jointStatistic.getLarger1Count(q+1);
        sumAXBhalf -= countsAXBhalf[q];
        sumBXAhalf -= countsBXAhalf[q];
    }
    countsAXBhalf[jointStatistic.getQ()] = sumAXBhalf;
    countsBXAhalf[jointStatistic.getQ()] = sumBXAhalf;
    const MaxLikelihoodEstimator estimator2(jointStatistic.getP(), jointStatistic.getQ()-1);

    const double cardinalityAXBhalf = estimator2(countsAXBhalf);
    const double cardinalityBXAhalf = estimator2(countsBXAhalf);

    cardinalityA = cardinalityABX - cardinalityBX;
    cardinalityB = cardinalityABX - cardinalityAX;
    double cardinalityX1 = 1.5*cardinalityBX + 1.5*cardinalityAX - cardinalityBXAhalf - cardinalityAXBhalf;
    double cardinalityX2 = 2.*(cardinalityBXAhalf + cardinalityAXBhalf) - 3*cardinalityABX;

    cardinalityX = std::max(0., 0.5*(cardinalityX1 + cardinalityX2));
}




void inclusionExclusionTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX) {

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    const double cardinalityAX = estimator(jointStatistic.get1Counts());
    const double cardinalityBX = estimator(jointStatistic.get2Counts());
    const double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    cardinalityA = cardinalityABX - cardinalityBX;
    cardinalityB = cardinalityABX - cardinalityAX;
    cardinalityX = std::max(0., cardinalityBX + cardinalityAX - cardinalityABX);
}


class JointLogLikelihoodFunction {

    const TwoHyperLogLogStatistic& jointStatistic;
    const int q;
    std::vector<double> pow2k;
    std::size_t* numFunctionEvaluations;
    std::size_t* numGradientEvaluations;

    double linTermA, linTermB, linTermX;

    static void calcExp(double x, double& y, double& z, double& xy) {
        static const double ln2 = std::log(2);
        if (x >= ln2) {
            y = std::exp(-x);
            z = 1 - y;
        }
        else {
            z = -std::expm1(-x);
            y = 1 - z;
        }
        if (y > 0) {
            xy = x * y;
        }
        else {
            xy = 0;
        }
    }

public:

    void evaluate(
        const double phiA,
        const double phiB,
        const double phiX,
        const bool calculateGradient,
        double& f,
        double& fa,
        double& fb,
        double& fx) const {

        (*numFunctionEvaluations) += 1;
        if (calculateGradient) {
            (*numGradientEvaluations) += 1;
        }

        const double expPhiA = std::exp(phiA);
        const double expPhiB = std::exp(phiB);
        const double expPhiX = std::exp(phiX);

        f = 0.;
        fa = 0.;
        fb = 0.;
        fx = 0.;

        for (int k = q + 1; k >= 1; --k) {

            int cSmaller1 = jointStatistic.getSmaller1Count(k);
            int cLarger1 = jointStatistic.getLarger1Count(k);
            int cSmaller2 = jointStatistic.getSmaller2Count(k);
            int cLarger2 = jointStatistic.getLarger2Count(k);
            int cEqual = jointStatistic.getEqualCount(k);

            double ya = 0, za = 0, xaya = 0;
            double yb = 0, zb = 0, xbyb = 0;
            double yx = 0, zx = 0, xxyx = 0;

            if (cSmaller1 > 0 || cEqual > 0 || cLarger1 > 0) {
                double xa = expPhiA * pow2k[k];
                calcExp(xa, ya, za, xaya);
            }
            if (cSmaller2 > 0 || cEqual > 0 || cLarger2 > 0) {
                double xb = expPhiB * pow2k[k];
                calcExp(xb, yb, zb, xbyb);
            }
            if (cSmaller1 > 0 || cEqual > 0 || cSmaller2 > 0) {
                double xx = expPhiX * pow2k[k];
                calcExp(xx, yx, zx, xxyx);
            }

            if (cSmaller1 > 0) {
                double arg = zx + yx * za;
                f  += cSmaller1 * std::log(arg);
                if (calculateGradient) {
                    double tmp = cSmaller1 / arg;
                    fa += tmp * yx * xaya;
                    fx += tmp * ya * xxyx;
                }
            }

            if (cLarger1 > 0) {
                f  += cLarger1 * std::log(za);
                if (calculateGradient) fa += cLarger1 * xaya / za;
            }

            if (cSmaller2 > 0) {
                double arg = zx + yx * zb;
                f  += cSmaller2 * std::log(arg);
                if (calculateGradient) {
                    double tmp = cSmaller2 / arg;
                    fb += tmp * yx * xbyb;
                    fx += tmp * yb * xxyx;
                }
            }

            if (cLarger2 > 0) {
                f  += cLarger2 * std::log(zb);
                if (calculateGradient) fb += cLarger2 * xbyb / zb;
            }

            if (cEqual > 0) {
                double arg = za * zb * yx + zx;
                f  += cEqual * std::log(arg);
                if (calculateGradient) {
                    double tmp = cEqual / arg;
                    double tmpyx = tmp * yx;
                    fa += tmpyx * zb * xaya;
                    fb += tmpyx * za * xbyb;
                    fx += tmp * (ya + yb * za) * xxyx;
                }
            }
        }

        double linTermExpPhiA = linTermA * expPhiA;
        double linTermExpPhiB = linTermB * expPhiB;
        double linTermExpPhiX = linTermX * expPhiX;

        f  -= linTermExpPhiA + linTermExpPhiB + linTermExpPhiX;
        fa -= linTermExpPhiA;
        fb -= linTermExpPhiB;
        fx -= linTermExpPhiX;
    }

public:

    JointLogLikelihoodFunction(const TwoHyperLogLogStatistic& jointStatistic_) :
        jointStatistic(jointStatistic_),
        q(jointStatistic_.getQ()),
        pow2k(q+2),
        numFunctionEvaluations(new std::size_t(0)),
        numGradientEvaluations(new std::size_t(0)) {

        pow2k[q] = ldexp(1., -q-jointStatistic.getP());
        pow2k[q+1] = pow2k[q];
        for (int k = q; k >= 1; --k) {
            pow2k[k-1] = pow2k[k] + pow2k[k];
        }

        linTermA = 0;
        linTermB = 0;
        linTermX = 0;
        for (int k = q; k >= 1; --k) {
            int cSmaller1 = jointStatistic.getSmaller1Count(k);
            int cLarger1 = jointStatistic.getLarger1Count(k);
            int cSmaller2 = jointStatistic.getSmaller2Count(k);
            int cLarger2 = jointStatistic.getLarger2Count(k);
            int cEqual = jointStatistic.getEqualCount(k);
            linTermA += cSmaller1 + cEqual + cLarger1;
            linTermB += cSmaller2 + cEqual + cLarger2;
            linTermX += cSmaller1 + cEqual + cSmaller2;
            linTermA *= 0.5;
            linTermB *= 0.5;
            linTermX *= 0.5;
        }
        linTermA += jointStatistic.get1Count(0);
        linTermB += jointStatistic.get2Count(0);
        linTermX += jointStatistic.getMinCount(0);
        linTermA /= jointStatistic.getNumRegisters();
        linTermB /= jointStatistic.getNumRegisters();
        linTermX /= jointStatistic.getNumRegisters();
    }

    ~JointLogLikelihoodFunction() {
        delete numFunctionEvaluations;
        delete numGradientEvaluations;
    }

    std::size_t getNumFunctionEvaluations() const {
        return *numFunctionEvaluations;
    }

    std::size_t getNumGradientEvaluations() const {
        return *numGradientEvaluations;
    }
};

void estimateInitialCardinalities(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, bool& isOptimum) {

    isOptimum = false;

    const int m = jointStatistic.getNumRegisters();

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());
    const double cardinalityAX = estimator(jointStatistic.get1Counts());
    const double cardinalityBX = estimator(jointStatistic.get2Counts());

    if(jointStatistic.getMinCount(0) == m) {
        cardinalityX = 0;
        cardinalityA = cardinalityAX;
        cardinalityB = cardinalityBX;
        isOptimum = true;
        return;
    }

    const std::vector<int>& larger1Counts = jointStatistic.getLarger1Counts();
    if(std::all_of(larger1Counts.begin(), larger1Counts.end(), [](int i) { return i==0; })) {
        cardinalityA = 0;
        cardinalityB = std::max(0., cardinalityBX - cardinalityAX);
        cardinalityX = cardinalityAX;
        return;
    }
    const std::vector<int>& larger2Counts = jointStatistic.getLarger2Counts();
    if(std::all_of(larger2Counts.begin(), larger2Counts.end(), [](int i) { return i==0; })) {
        cardinalityB = 0;
        cardinalityA = std::max(0., cardinalityAX - cardinalityBX);
        cardinalityX = cardinalityBX;
        return;
    }

    const double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    cardinalityA = std::max(1., cardinalityABX - cardinalityBX);
    cardinalityB = std::max(1., cardinalityABX - cardinalityAX);
    cardinalityX = std::max(1., cardinalityAX + cardinalityBX - cardinalityABX);

}

#ifdef CARDINALITY_ESTIMATION_USE_CERES

class LogLikelihoodFunctionForCeres : public ceres::FirstOrderFunction {

    const JointLogLikelihoodFunction& _logLikelihoodFunction;

public:

    LogLikelihoodFunctionForCeres(const JointLogLikelihoodFunction& logLikelihoodFunction) : _logLikelihoodFunction(logLikelihoodFunction) {}

    virtual bool Evaluate(const double* parameters, double* cost, double* gradient) const {

        const double phiA = parameters[0];
        const double phiB = parameters[1];
        const double phiX = parameters[2];

        double f, fa, fb, fx;
        _logLikelihoodFunction.evaluate(phiA, phiB, phiX, gradient != NULL, f, fa, fb, fx);

        cost[0] = -f;
        if (gradient != NULL) {
            gradient[0] = -fa;
            gradient[1] = -fb;
            gradient[2] = -fx;
        }
        return true;
    }

    virtual int NumParameters() const { return 3; }

};

class TerminationCheckingCallback : public ceres::IterationCallback {
    const double eps_;
public:
    TerminationCheckingCallback(double eps ) : eps_(eps) {}

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) {

        if (summary.step_is_successful && summary.step_norm <= eps_) {
            return ceres::SOLVER_TERMINATE_SUCCESSFULLY;
        }
        return ceres::SOLVER_CONTINUE;
    }
};

void maxLikelihoodTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, std::size_t& numFunctionEvaluations, std::size_t& numGradientEvaluations) {

    const double eps = 1e-2;
    const int m = jointStatistic.getNumRegisters();
    const double relativeErrorLimit = eps/(sqrt(m));

    const double zeroPhi = 2 * std::log(std::numeric_limits<double>::min());
    assert(std::exp(zeroPhi) == 0);

    bool isOptimum;
    estimateInitialCardinalities(jointStatistic, cardinalityA, cardinalityB, cardinalityX, isOptimum);
    if (isOptimum) {
        numFunctionEvaluations = 0;
        numGradientEvaluations = 0;
        return;
    }

    double phiA = (cardinalityA > 0)?std::log(cardinalityA):zeroPhi;
    double phiB = (cardinalityB > 0)?std::log(cardinalityB):zeroPhi;
    double phiX = (cardinalityX > 0)?std::log(cardinalityX):zeroPhi;
    double parameters[3] = {phiA, phiB, phiX};

    JointLogLikelihoodFunction logLikelihoodFunction(jointStatistic);

    ceres::GradientProblem problem(new LogLikelihoodFunctionForCeres(logLikelihoodFunction));
    ceres::GradientProblemSolver::Options options;
    options.line_search_direction_type = ceres::BFGS;
    options.function_tolerance = 0;
    options.gradient_tolerance = 0;
    options.parameter_tolerance = 0;
    options.logging_type = ceres::SILENT;
    TerminationCheckingCallback callback(relativeErrorLimit);
    options.callbacks.push_back(&callback);

    ceres::GradientProblemSolver::Summary summary;
    ceres::Solve(options, problem, parameters, &summary);
    numFunctionEvaluations = logLikelihoodFunction.getNumFunctionEvaluations();
    numGradientEvaluations = logLikelihoodFunction.getNumGradientEvaluations();
    assert(summary.termination_type == ceres::USER_SUCCESS || summary.termination_type == ceres::CONVERGENCE);

    cardinalityA = std::exp(parameters[0]);
    cardinalityB = std::exp(parameters[1]);
    cardinalityX = std::exp(parameters[2]);
}

#endif

#ifdef CARDINALITY_ESTIMATION_USE_DLIB


class LogLikelihoodFunctionForDlib {

    const JointLogLikelihoodFunction& _logLikelihoodFunction;
    double fa;
    double fb;
    double fx;
    bool status = false;

public:
    LogLikelihoodFunctionForDlib(const JointLogLikelihoodFunction& logLikelihoodFunction) :
        _logLikelihoodFunction(logLikelihoodFunction) {}

    double value(const dlib::matrix<double,3,1>& x) {
        assert(status == false); // check if value and derivative are alternately called
        double f;
        _logLikelihoodFunction.evaluate(x(0), x(1), x(2), true, f, fa, fb, fx);
        status = true;
        return f;
    }

    dlib::matrix<double,3,1> derivative(const dlib::matrix<double,3,1>& x) {
        assert(status == true); // check if value and derivative are alternately called
        status = false;
        return {fa, fb, fx};
    }

};

class StopStrategy {
    const double _min_delta;
    bool _been_used;
    dlib::matrix<double,3,1> _previous_values;
public:

    StopStrategy(double min_delta) : _min_delta(min_delta),  _been_used(false) {}

    bool should_continue_search (const dlib::matrix<double,3,1>& values, const double, const dlib::matrix<double,3,1>&)
    {
        bool ret =
            !_been_used ||
            std::abs(values(0) - _previous_values(0)) > _min_delta ||
            std::abs(values(1) - _previous_values(1)) > _min_delta ||
            std::abs(values(2) - _previous_values(2)) > _min_delta;
        _been_used = true;
        _previous_values = values;
        return ret;
    }
};


void maxLikelihoodTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, std::size_t& numFunctionEvaluations, std::size_t& numGradientEvaluations) {

    const double eps = 1e-2;

    const int m = jointStatistic.getNumRegisters();
    const double relativeErrorLimit = eps/(sqrt(m));

    const double zeroPhi = 2 * std::log(std::numeric_limits<double>::min());
    assert(std::exp(zeroPhi) == 0);

    bool isOptimum;
    estimateInitialCardinalities(jointStatistic, cardinalityA, cardinalityB, cardinalityX, isOptimum);
    if (isOptimum) {
        numFunctionEvaluations = 0;
        numGradientEvaluations = 0;
        return;
    }

    double phiA = (cardinalityA > 0)?std::log(cardinalityA):zeroPhi;
    double phiB = (cardinalityB > 0)?std::log(cardinalityB):zeroPhi;
    double phiX = (cardinalityX > 0)?std::log(cardinalityX):zeroPhi;

    JointLogLikelihoodFunction logLikelihoodFunction(jointStatistic);

    LogLikelihoodFunctionForDlib logLikelihoodFunctionForDlib(logLikelihoodFunction);

    dlib::matrix<double,3,1> starting_point = {phiA, phiB, phiX};
    dlib::find_max(
        dlib::bfgs_search_strategy(),  // Use BFGS search algorithm
        StopStrategy(relativeErrorLimit),
        std::bind(&LogLikelihoodFunctionForDlib::value, std::ref(logLikelihoodFunctionForDlib), std::placeholders::_1),
        std::bind(&LogLikelihoodFunctionForDlib::derivative, std::ref(logLikelihoodFunctionForDlib), std::placeholders::_1),
        starting_point,
        std::numeric_limits<double>::infinity());

    cardinalityA = std::exp(starting_point(0));
    cardinalityB = std::exp(starting_point(1));
    cardinalityX = std::exp(starting_point(2));
    numFunctionEvaluations = logLikelihoodFunction.getNumFunctionEvaluations();
    numGradientEvaluations = logLikelihoodFunction.getNumGradientEvaluations();
}

#endif

#endif // _CARDINALITY_ESTIMATION_HPP_
