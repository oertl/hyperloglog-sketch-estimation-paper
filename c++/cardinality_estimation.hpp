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

class OptimalLinearCountingEstimator {

    const int p;
    const int m = 1 << p;
    const std::vector<double> results = initResults(m);

    static std::vector<double> initResults(int m) {

        std::vector<double> results(m+1);
        results[m] = 0;
        for (int c = m-1; c >= 0; --c) {
            results[c] = results[c + 1] + m/(double)c;
        }
        return results;
    }

public:

    OptimalLinearCountingEstimator(const int p_) : p(p_) {}

    double operator()(const std::vector<int>& c) const {
        return results[c[0]];
    }
};


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

// the minimum of this function gives the max likelihood estimate
void eval_joint_log_likelihood_function_and_derivatives(
    const TwoHyperLogLogStatistic& jointStatistic,
    const double phiA,
    const double phiB,
    const double phiX,
    double& f,
    double& fa,
    double& fb,
    double& fx) {

    const int q = jointStatistic.getQ();

    const double expPhiA = std::exp(phiA);
    const double expPhiB = std::exp(phiB);
    const double expPhiX = std::exp(phiX);

    const double expPhiAdivAX = expPhiA / (expPhiA + expPhiX);   // TODO handle division by 0?
    const double expPhiXdivAX = expPhiX / (expPhiA + expPhiX);
    const double expPhiBdivBX = expPhiB / (expPhiB + expPhiX);
    const double expPhiXdivBX = expPhiX / (expPhiB + expPhiX);


    double expm1a[q+1];  // 1-exp(-exp(phiA/2^k))
    double expm1b[q+1];  // 1-exp(-exp(phiB/2^k))
    double expm1ax[q+1]; // 1-exp(-exp(phiA/2^k)-exp(phiX/2^k))
    double expm1bx[q+1]; // 1-exp(-exp(phiB/2^k)-exp(phiX/2^k))
    double expa[q+1];  // exp(-exp(phiA/2^k))
    double expb[q+1];  // exp(-exp(phiB/2^k))
    double expax[q+1]; // exp(-exp(phiA/2^k)-exp(phiX/2^k))
    double expbx[q+1]; // exp(-exp(phiB/2^k)-exp(phiX/2^k))
    double expm1x[q+1];
    double expx[q+1];
    double expPhiApow2k[q+1]; // exp(phiA/2^k)
    double expPhiBpow2k[q+1]; // exp(phiB/2^k)
    double expPhiXpow2k[q+1]; // exp(phiX/2^k)
    double expPhiAXpow2k[q+1]; //exp(phiA/2^k) + exp(phiX/2^k)
    double expPhiBXpow2k[q+1]; //exp(phiB/2^k) + exp(phiX/2^k)

    expPhiApow2k[0] = expPhiA;
    expPhiBpow2k[0] = expPhiB;
    expPhiXpow2k[0] = expPhiX;
    for (int k = 1; k <= q; ++k) {
        expPhiApow2k[k] = expPhiApow2k[k-1] * 0.5;
        expPhiBpow2k[k] = expPhiBpow2k[k-1] * 0.5;
        expPhiXpow2k[k] = expPhiXpow2k[k-1] * 0.5;
    }

    for (int k = 0; k <= q; ++k) {
        expPhiAXpow2k[k] = expPhiApow2k[k] + expPhiXpow2k[k];
        expPhiBXpow2k[k] = expPhiBpow2k[k] + expPhiXpow2k[k];
    }


    for (int k = q; k >= 1; --k) {

        const int cEqual    = jointStatistic.getEqualCounts()[k];
        const int cLarger2  = jointStatistic.getLarger2Counts()[k];
        const int cSmaller1 = jointStatistic.getSmaller1Counts()[k];
        const int cSmaller2 = jointStatistic.getSmaller2Counts()[k];
        const int cLarger1  = jointStatistic.getLarger1Counts()[k];

        if (cLarger1 > 0 || cEqual > 0) {
            expm1a[k] = -std::expm1(-expPhiApow2k[k]);
            expa[k] = std::exp(-expPhiApow2k[k]);
        }
        if (cLarger2 > 0 || cEqual > 0) {
            expm1b[k] = -std::expm1(-expPhiBpow2k[k]);
            expb[k] = std::exp(-expPhiBpow2k[k]);
        }
        if (cSmaller1 > 0 || cEqual > 0) {
            expm1ax[k] = -std::expm1(-expPhiAXpow2k[k]);
            expax[k] = std::exp(-expPhiAXpow2k[k]);
        }
        if (cSmaller2 > 0 || cEqual > 0) {
            expm1bx[k] = -std::expm1(-expPhiBXpow2k[k]);
            expbx[k] = std::exp(-expPhiBXpow2k[k]);
        }

        if (cEqual > 0) {
            expm1x[k] = -std::expm1(-expPhiXpow2k[k]);
            expx[k] = std::exp(-expPhiXpow2k[k]);
        }
    }

    f = 0.;
    fa = 0.;
    fb = 0.;
    fx = 0.;

    double term1a = 0;
    double term1b = 0;
    double term1x = 0;

    for (int k = q; k >= 0; --k) {
        const int cEqual    = jointStatistic.getEqualCounts()[k];
        const int cLarger2  = jointStatistic.getLarger2Counts()[k];
        const int cSmaller1 = jointStatistic.getSmaller1Counts()[k];
        const int cSmaller2 = jointStatistic.getSmaller2Counts()[k];
        const int cLarger1  = jointStatistic.getLarger1Counts()[k];
        term1a += (cEqual + cSmaller1 + cLarger1 ) * expPhiApow2k[k];
        term1b += (cEqual + cLarger2  + cSmaller2) * expPhiBpow2k[k];
        term1x += (cEqual + cSmaller1 + cSmaller2) * expPhiXpow2k[k];
    }

    for (int kk = q+1; kk >= 1; --kk) {
        const int k = std::min(kk, q);
        const int cEqual    = jointStatistic.getEqualCounts()[k];
        const int cLarger2  = jointStatistic.getLarger2Counts()[k];
        const int cSmaller1 = jointStatistic.getSmaller1Counts()[k];
        const int cSmaller2 = jointStatistic.getSmaller2Counts()[k];
        const int cLarger1  = jointStatistic.getLarger1Counts()[k];

        if (cSmaller1 > 0) {
            f  -= cSmaller1 * std::log(expm1ax[k]);
            double tmp = cSmaller1 * expax[k] * (expPhiAXpow2k[k] / expm1ax[k]);
            fa -= tmp * expPhiAdivAX;
            fx -= tmp * expPhiXdivAX;
        }

        if (cLarger1 > 0) {
            f  -= cLarger1 * std::log(expm1a[k]);
            fa -= cLarger1 * expa[k] * (expPhiApow2k[k] / expm1a[k]);
        }

        if (cSmaller2 > 0) {
            f  -= cSmaller2 * std::log(expm1bx[k]);
            double tmp = cSmaller2 * expbx[k] * (expPhiBXpow2k[k] / expm1bx[k]);
            fb -= tmp * expPhiBdivBX;
            fx -= tmp * expPhiXdivBX;
        }

        if (cLarger2 > 0) {
            f  -= cLarger2 * std::log(expm1b[k]);
            fb -= cLarger2 * expb[k] * (expPhiBpow2k[k] / expm1b[k]);
        }

        if (cEqual > 0) {
            f  -= cEqual * std::log(expm1a[k] * expm1b[k] * expx[k] + expm1x[k]);
            fa -= cEqual * (((expax[k] * expm1b[k])/(expm1a[k] * expm1b[k] * expx[k] + expm1x[k])) * expPhiApow2k[k]);
            fb -= cEqual * (((expbx[k] * expm1a[k])/(expm1a[k] * expm1b[k] * expx[k] + expm1x[k])) * expPhiBpow2k[k]);
            fx -= cEqual * (((expx[k] *(1 - expm1a[k] * expm1b[k]))/(expm1a[k] * expm1b[k] * expx[k] + expm1x[k])) * expPhiXpow2k[k]);
        }

    }

    f += term1a + term1b + term1x;
    fa += term1a;
    fb += term1b;
    fx += term1x;

}

void log_likelihood_function_value_and_derivatives_for_gsl(const gsl_vector *phi, void *params, double *f, gsl_vector *fGrad)
{

    const double phiA = gsl_vector_get(phi, 0);
    const double phiB = gsl_vector_get(phi, 1);
    const double phiX = gsl_vector_get(phi, 2);

    const TwoHyperLogLogStatistic& jointStatistic = *((TwoHyperLogLogStatistic*)params);

    double fa, fb, fx;

    eval_joint_log_likelihood_function_and_derivatives(
        jointStatistic,
        phiA, phiB, phiX,
        *f, fa, fb, fx);

    gsl_vector_set(fGrad, 0, fa);
    gsl_vector_set(fGrad, 1, fb);
    gsl_vector_set(fGrad, 2, fx);

}

void log_likelihood_function_derivatives_for_gsl(const gsl_vector *phi, void *params, gsl_vector *fGrad)
{
    double f;
    log_likelihood_function_value_and_derivatives_for_gsl(phi, params, &f, fGrad);
}

double log_likelihood_function_value_for_gsl(const gsl_vector *phi, void *params)
{
    gsl_vector *fGrad = gsl_vector_alloc(3);
    double f;
    log_likelihood_function_value_and_derivatives_for_gsl(phi, params, &f, fGrad);
    gsl_vector_free (fGrad);
    return f;
}

void maxLikelihoodTwoHyperLogLogEstimation(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, bool& maxNumIterationsReached, bool& iterationAborted, int& numIterations) {

    const double eps = 1e-5;

    const int maxNumIterations = 100000;

    maxNumIterationsReached = false;
    iterationAborted = false;

    const int m = jointStatistic.getNumRegisters();

    // set initial vector
    // TODO improve initial values

    // TODO add handling of special cases (all register in 1 are smaller than in 2, etc.)

    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    //const double cardinalityAX = estimator(jointStatistic.get1Counts());
    //const double cardinalityBX = estimator(jointStatistic.get2Counts());
    //const double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    //const double cardinalityMin = estimator(jointStatistic.getMinCounts());

    const double initalStepFactor = 2.;

    double initCadinalityA = 1.;
    double initCadinalityB = 1.;
    double initCadinalityX = 1.;

    double lastPhiA = std::log(initCadinalityA/m);
    double lastPhiB = std::log(initCadinalityB/m);
    double lastPhiX = std::log(initCadinalityX/m);
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
            std::fabs(currentPhiA - lastPhiA) <= eps &&
            std::fabs(currentPhiB - lastPhiB) <= eps &&
            std::fabs(currentPhiX - lastPhiX) <= eps) {
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

void analyticalJointHyperLogLogEstimator(const TwoHyperLogLogStatistic& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX) {

    int q = jointStatistic.getQ();
    int m = jointStatistic.getNumRegisters();

    const std::vector<int>& pmfAX = jointStatistic.get1Counts();
    const std::vector<int>& pmfBX = jointStatistic.get2Counts();
    const std::vector<int>& pmfABX = jointStatistic.getMaxCounts();

    double cdfAX[q+1];
    double cdfBX[q+1];
    double cdfABX[q+1];

    int sumAX = 0;
    int sumBX = 0;
    int sumABX = 0;
    for(int i = 0; i < q+1; ++i) {
        sumAX += pmfAX[i];
        cdfAX[i] = sumAX/(double)m;
        sumBX += pmfBX[i];
        cdfBX[i] = sumBX/(double)m;
        sumABX += pmfABX[i];
        cdfABX[i] = sumABX/(double)m;
     }

    // estimate cdfs of A, B, X
    std::vector<double> cdfA(q+1);
    std::vector<double> cdfB(q+1);
    std::vector<double> cdfX(q+1);
    for(int i = q; i >=0; --i) {
        if (cdfABX[i] > 0) {
            cdfX[i] = cdfAX[i] * cdfBX[i] / cdfABX[i];
            assert(cdfX[i] <= 2);
        }
        else {
            if (i + 1 <= q) {
                cdfX[i] = cdfX[i+1] * cdfX[i+1];
            }
            else {
                assert(false);
                cdfX[i] = 1; // TODO
            }
        }
    }
    for(int i = q; i >=0; --i) {
        if (cdfBX[i] > 0) {
            cdfA[i] = cdfABX[i]/cdfBX[i];
            assert(cdfA[i] <= 2);
        }
        else {
            if (i + 1 <= q) {
                cdfA[i] = cdfA[i+1] * cdfA[i+1];
            }
            else {
                assert(false);
                cdfA[i] = 1; // TODO
            }
        }
    }
    for(int i = q; i >=0; --i) {
        if (cdfAX[i] > 0) {
            cdfB[i] = cdfABX[i]/cdfAX[i];
            assert(cdfB[i] <= 2);
        }
        else {
            if (i + 1 <= q) {
                cdfB[i] = cdfB[i+1] * cdfB[i+1];
            }
            else {
                assert(false);
                cdfB[i] = 1; // TODO
            }
        }
    }

    cardinalityA = estimateCardinalityFromCdf(cdfA, m);
    cardinalityB = estimateCardinalityFromCdf(cdfB, m);
    cardinalityX = estimateCardinalityFromCdf(cdfX, m);
}


#endif // _CARDINALITY_ESTIMATION_HPP_
