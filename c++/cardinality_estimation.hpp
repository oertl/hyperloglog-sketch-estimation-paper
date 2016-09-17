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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

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


double estimateCardinalityFromCdf(const std::vector<double> cdf, int m) {
    int q = cdf.size() - 1;

    if (cdf[0] >= 1.) std::numeric_limits<double>::infinity();

    double lowerTail = 0;
    {
        double arg = cdf[0];
        double pow2 = 1;
        while(true) {
            const double lowerTailBak = lowerTail;
            lowerTail += pow2 * arg;
            if (lowerTail == lowerTailBak) break;
            pow2 += pow2;
            arg *= arg;
        }
    }

    double upperTail = 0;
    {
        double arg = std::sqrt(cdf[cdf.size()-1]);
        double pow2 = 0.5;
        while(true) {
            const double upperTailBak = upperTail;
            upperTail += pow2 * arg;
            if (upperTail == upperTailBak) break;
            pow2 *= 0.5;
            arg = std::sqrt(arg);
        }
    }

    double sum = upperTail;
    for (int k = q; k >= 1; --k) {
        sum += cdf[k];
        sum *= 0.5;
    }

    sum += lowerTail;

    return (m/std::log(2))/sum;

}

class MaxLikelihoodEstimator {
    const double eps = 1e-2;
    const int p;
    const int q;
    const int m;
    const double relativeErrorLimit;

public:

    MaxLikelihoodEstimator(const int p_, const int q_) :p(p_), q(q_), m(1 << p_), relativeErrorLimit(eps/(sqrt(m))) {}

    double estimate(const std::vector<int>& c, int& outerLoopIterationsCount, int& innerLoop1IterationsCount, int& innerLoop2IterationsCount, int& logEvaluationCount, int& kMin, int& kMax) const {

        outerLoopIterationsCount = 0;
        innerLoop1IterationsCount = 0;
        innerLoop2IterationsCount = 0;
        logEvaluationCount = 0;

        for(kMin=0; c[kMin]==0; ++kMin);
        int kMinPrime = std::max(1, kMin);

        for(kMax=q+1; c[kMax]==0; --kMax);
        int kMaxPrime = std::min(q, kMax);

        if (kMin > q) return std::numeric_limits<double>::infinity();

        double z = 0.;
        int mPrime = c[q+1];
        double y = ldexp(1., -kMaxPrime);
        for (int k = kMaxPrime; k >= kMinPrime; --k) {
            z += y*c[k];
            y += y;
            mPrime += c[k];
        }

        int cPrime = c[q+1];
        if (q >= 1) {
            cPrime += c[kMaxPrime];
        }

        double gPrevious = 0;

        double x;
        double a = z + c[0];
        {
            double b = z + ldexp(c[q+1], -q);
            if (b <= 1.5*a) { // TODO constant 1.5
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
        int kMin;
        int kMax;

        return estimate(c, outerLoopIterationsCount, innerLoop1IterationsCount, innerLoop2IterationsCount, logEvaluationCount, kMin, kMax);
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

class CorrectedRawEstimator {
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
        } while(zPrime < z);
        return z;
    }

    static double tau(double x, int& numIterations) {
        numIterations = 0;
        double zPrime;
        double y = 1.0;
        double z = 0;
        do {
            numIterations += 1;
            x = std::sqrt(x);
            zPrime = z;
            y *= 0.5;
            z += (1 - x)*x*y;
        } while(zPrime < z);
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

    CorrectedRawEstimator(const int p_, const int q_) : p(p_), q(q_) {}

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

void inclusionExclusionTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX) {

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    const double cardinalityAX = estimator(jointStatistic.get1Counts());
    const double cardinalityBX = estimator(jointStatistic.get2Counts());
    const double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    cardinalityA = cardinalityABX - cardinalityBX;
    cardinalityB = cardinalityABX - cardinalityAX;
    cardinalityX = std::max(0., cardinalityBX + cardinalityAX - cardinalityABX);
}

void calcExp(double x, double& a, double& b) {
    static const double ln2 = std::log(2);
    if (x >= ln2) {
        a = std::exp(-x);
        b = 1 - a;
    }
    else {
        b = -std::expm1(-x);
        a = 1 - b;
    }
}

// the minimum of this function gives the max likelihood estimate
void eval_joint_log_likelihood_function_and_derivatives(
    const TwoHyperLogLogStatistic& jointStatistic,
    const double phiA,
    const double phiB,
    const double phiX,
    double& f,
    double& fa,
    double& fb,
    double& fx,
    bool calcValue,
    bool calcDerivative) {

    const int q = jointStatistic.getQ();

    const double expPhiA = std::exp(phiA);
    const double expPhiB = std::exp(phiB);
    const double expPhiX = std::exp(phiX);

    f = 0.;
    fa = 0.;
    fb = 0.;
    fx = 0.;

    double linTermA = 0;
    double linTermB = 0;
    double linTermX = 0;

    double pow2k = ldexp(1., -q);
    for (int k = q + 1; k >= 1; --k) {

        int cSmaller1 = jointStatistic.getSmaller1Count(k);
        int cLarger1 = jointStatistic.getLarger1Count(k);
        int cSmaller2 = jointStatistic.getSmaller2Count(k);
        int cLarger2 = jointStatistic.getLarger2Count(k);
        int cEqual = jointStatistic.getEqualCount(k);

        double xak, yak, zak;
        double xbk, ybk, zbk;
        double xxk, yxk, zxk;

        if (cSmaller1 > 0 || cEqual > 0 || cLarger1 > 0) {
            xak = expPhiA * pow2k;
            calcExp(xak, yak, zak);
        }
        if (cSmaller2 > 0 || cEqual > 0 || cLarger2 > 0) {
            xbk = expPhiB * pow2k;
            calcExp(xbk, ybk, zbk);
        }
        if (cSmaller1 > 0 || cEqual > 0 || cSmaller2 > 0) {
            xxk = expPhiX * pow2k;
            calcExp(xxk, yxk, zxk);
        }

        if (cSmaller1 > 0) {
            double arg = zxk + yxk*zak;
            if (calcValue) f  -= cSmaller1 * std::log(arg);
            if (calcDerivative) {
                double tmp = cSmaller1 * yak * yxk / arg;
                fa -= tmp * xak;
                fx -= tmp * xxk;
            }
        }

        if (cLarger1 > 0) {
            if (calcValue) f  -= cLarger1 * std::log(zak);
            if (calcDerivative) fa -= cLarger1 * yak * xak / zak;
        }

        if (cSmaller2 > 0) {
            double arg = zxk + yxk*zbk;
            if (calcValue) f  -= cSmaller2 * std::log(arg);
            if (calcDerivative) {
                double tmp = cSmaller2 * ybk * yxk / arg;
                fb -= tmp * xbk;
                fx -= tmp * xxk;
            }
        }

        if (cLarger2 > 0) {
            if (calcValue) f  -= cLarger2 * std::log(zbk);
            if (calcDerivative) fb -= cLarger2 * ybk * xbk / zbk;
        }

        if (cEqual > 0) {
            double arg = zak * zbk * yxk + zxk;
            if (calcValue) f  -= cEqual * std::log(arg);
            if (calcDerivative) {
                double yazb = yak*zbk;
                double ybza = ybk*zak;
                double tmp = cEqual * yxk / arg;
                fa -= tmp * yazb * xak;
                fb -= tmp * ybza * xbk;
                fx -= tmp * (yazb + ybza + yak * ybk) * xxk;
            }
        }
        if (k <= q) {
            linTermA += cSmaller1 + cEqual + cLarger1;
            linTermB += cSmaller2 + cEqual + cLarger2;
            linTermX += cSmaller1 + cEqual + cSmaller2;
            linTermA *= 0.5;
            linTermB *= 0.5;
            linTermX *= 0.5;
            pow2k += pow2k;
        }
    }

    linTermA += jointStatistic.get1Count(0);
    linTermB += jointStatistic.get2Count(0);
    linTermX += jointStatistic.getMinCount(0);
    linTermA *= expPhiA;
    linTermB *= expPhiB;
    linTermX *= expPhiX;

    f += linTermA + linTermB + linTermX;
    fa += linTermA;
    fb += linTermB;
    fx += linTermX;

}

void log_likelihood_function_value_and_derivatives_for_gsl(const gsl_vector *phi, void *params, double *f, gsl_vector *fGrad)
{

    const TwoHyperLogLogStatistic& jointStatistic = *((TwoHyperLogLogStatistic*)params);

    double fa, fb, fx;

    eval_joint_log_likelihood_function_and_derivatives(
        jointStatistic,
        gsl_vector_get(phi, 0), gsl_vector_get(phi, 1), gsl_vector_get(phi, 2),
        *f, fa, fb, fx, true, true);

    gsl_vector_set(fGrad, 0, fa);
    gsl_vector_set(fGrad, 1, fb);
    gsl_vector_set(fGrad, 2, fx);

}

void log_likelihood_function_derivatives_for_gsl(const gsl_vector *phi, void *params, gsl_vector *fGrad)
{
    const TwoHyperLogLogStatistic& jointStatistic = *((TwoHyperLogLogStatistic*)params);

    double fa, fb, fx, f;

    eval_joint_log_likelihood_function_and_derivatives(
        jointStatistic,
        gsl_vector_get(phi, 0), gsl_vector_get(phi, 1), gsl_vector_get(phi, 2),
        f, fa, fb, fx, false, true);

    gsl_vector_set(fGrad, 0, fa);
    gsl_vector_set(fGrad, 1, fb);
    gsl_vector_set(fGrad, 2, fx);
}

double log_likelihood_function_value_for_gsl(const gsl_vector *phi, void *params)
{
    const TwoHyperLogLogStatistic& jointStatistic = *((TwoHyperLogLogStatistic*)params);

    double fa, fb, fx, f;

    eval_joint_log_likelihood_function_and_derivatives(
        jointStatistic,
        gsl_vector_get(phi, 0), gsl_vector_get(phi, 1), gsl_vector_get(phi, 2),
        f, fa, fb, fx, true, false);

    return f;
}

void maxLikelihoodTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, bool& maxNumIterationsReached, bool& iterationAborted, int& numIterations) {

    const double eps = 1e-2;
    const double initalStepFactor = 2;
    const int maxNumIterations = 10000;

    maxNumIterationsReached = false;
    iterationAborted = false;

    const int m = jointStatistic.getNumRegisters();
    const double relativeErrorLimit = eps/(sqrt(m));

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    double cardinalityAX = estimator(jointStatistic.get1Counts());
    double cardinalityBX = estimator(jointStatistic.get2Counts());

    // special handling, the sets A u X and B u X are disjoint, therefore X = O
    if (jointStatistic.getMinCount(0) == m) {
        cardinalityX = 0;
        cardinalityA = cardinalityAX;
        cardinalityB = cardinalityBX;
        return;
    }
    // for all remaining cases the maximum likelihood estimates for A,B,X will be positive

    double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    // set initial vector using inclusion exclusion principle
    double initCardinalityA = std::max(1., cardinalityABX - cardinalityBX);
    double initCardinalityB = std::max(1., cardinalityABX - cardinalityAX);
    double initCardinalityX = std::max(1., cardinalityAX + cardinalityBX - cardinalityABX);

    double lastPhiA = std::log(initCardinalityA/m);
    double lastPhiB = std::log(initCardinalityB/m);
    double lastPhiX = std::log(initCardinalityX/m);
    gsl_vector *phi = gsl_vector_alloc(3);
    gsl_vector_set(phi, 0, lastPhiA);
    gsl_vector_set(phi, 1, lastPhiB);
    gsl_vector_set(phi, 2, lastPhiX);

    gsl_multimin_function_fdf my_func;

    my_func.n = 3;
    my_func.f = log_likelihood_function_value_for_gsl;
    my_func.df = log_likelihood_function_derivatives_for_gsl;
    my_func.fdf = log_likelihood_function_value_and_derivatives_for_gsl;
    my_func.params = &const_cast<TwoHyperLogLogStatistic&>(jointStatistic);

    gsl_multimin_fdfminimizer *solver = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, 3); // 3 dimensions
    gsl_multimin_fdfminimizer_set(solver, &my_func, phi, std::log(initalStepFactor), 0.1);

    numIterations = 0;
    do
    {
        numIterations++;
        int status = gsl_multimin_fdfminimizer_iterate(solver);

        if (status) {
            if (status == GSL_ENOPROG) {
                iterationAborted = true;
                break;
            }
            std::cout << "error!" << std::endl;
            exit(-1);
        }

        const gsl_vector* currentPhi = gsl_multimin_fdfminimizer_x(solver);
        double currentPhiA = gsl_vector_get(currentPhi, 0);
        double currentPhiB = gsl_vector_get(currentPhi, 1);
        double currentPhiX = gsl_vector_get(currentPhi, 2);

        if (
            std::fabs(currentPhiA - lastPhiA) <= relativeErrorLimit &&
            std::fabs(currentPhiB - lastPhiB) <= relativeErrorLimit &&
            std::fabs(currentPhiX - lastPhiX) <= relativeErrorLimit) {
            break;
        }

        lastPhiA = currentPhiA;
        lastPhiB = currentPhiB;
        lastPhiX = currentPhiX;

        if(numIterations >= maxNumIterations) {
            maxNumIterationsReached = true;
            break;
        }
    }
    while (true);

    cardinalityA = std::exp(gsl_vector_get(solver->x, 0)) * m;
    cardinalityB = std::exp(gsl_vector_get(solver->x, 1)) * m;
    cardinalityX = std::exp(gsl_vector_get(solver->x, 2)) * m;

    gsl_multimin_fdfminimizer_free(solver);
    gsl_vector_free(phi);

}

#endif // _CARDINALITY_ESTIMATION_HPP_
