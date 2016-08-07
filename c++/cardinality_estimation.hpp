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

    f = 0.;
    fa = 0.;
    fb = 0.;
    fx = 0.;

    double term1a = 0;
    double term1b = 0;
    double term1x = 0;
    const double pow2q = std::ldexp(1., -q);
    double expPhiApow2k = expPhiA * pow2q;
    double expPhiBpow2k = expPhiB * pow2q;
    double expPhiXpow2k = expPhiX * pow2q;

    const double expPhiAdivAX = expPhiA / (expPhiA + expPhiX);   // TODO handle division by 0
    const double expPhiXdivAX = expPhiX / (expPhiA + expPhiX);
    const double expPhiBdivBX = expPhiB / (expPhiB + expPhiX);
    const double expPhiXdivBX = expPhiX / (expPhiB + expPhiX);

    for (int k = q; k >= 0; --k) {
        const int cEqual    = jointStatistic.getEqualCounts()[k];
        const int cLarger2  = jointStatistic.getLarger2Counts()[k];
        const int cSmaller1 = jointStatistic.getSmaller1Counts()[k];
        const int cSmaller2 = jointStatistic.getSmaller2Counts()[k];
        const int cLarger1  = jointStatistic.getLarger1Counts()[k];
        term1a += (cEqual + cSmaller1 + cLarger1 ) * expPhiApow2k;
        term1b += (cEqual + cLarger2  + cSmaller2) * expPhiBpow2k;
        term1x += (cEqual + cSmaller1 + cSmaller2) * expPhiXpow2k;
        expPhiApow2k += expPhiApow2k;
        expPhiBpow2k += expPhiBpow2k;
        expPhiXpow2k += expPhiXpow2k;
    }

    expPhiApow2k = expPhiA * pow2q;
    expPhiBpow2k = expPhiB * pow2q;
    expPhiXpow2k = expPhiX * pow2q;
    for (int k = q+1; k >= 1; --k) {
        const int cEqual    = jointStatistic.getEqualCounts()[k];
        const int cLarger2  = jointStatistic.getLarger2Counts()[k];
        const int cSmaller1 = jointStatistic.getSmaller1Counts()[k];
        const int cSmaller2 = jointStatistic.getSmaller2Counts()[k];
        const int cLarger1  = jointStatistic.getLarger1Counts()[k];

        double expm1a;
        double expm1b;
        double expm1ax;
        double expm1bx;

        const double expPhiAXpow2k = expPhiApow2k + expPhiXpow2k;
        const double expPhiBXpow2k = expPhiBpow2k + expPhiXpow2k;

        if (cLarger1 > 0 || cEqual > 0) {
            expm1a = -std::expm1(-expPhiApow2k);
        }
        if (cLarger2 > 0 || cEqual > 0) {
            expm1b = -std::expm1(-expPhiBpow2k);
        }
        if (cSmaller1 > 0 || cEqual > 0) {
            expm1ax = -std::expm1(-expPhiAXpow2k);
        }
        if (cSmaller2 > 0 || cEqual > 0) {
            expm1bx = -std::expm1(-expPhiBXpow2k);
        }

        if (cSmaller1 > 0) {
            f  -= cSmaller1 * std::log(expm1ax);
            double tmp = cSmaller1 * (1-expm1ax) * (expPhiAXpow2k / expm1ax);
            fa -= tmp * expPhiAdivAX;
            fx -= tmp * expPhiXdivAX;
        }

        if (cLarger1 > 0) {
            f  -= cLarger1 * std::log(expm1a);
            fa -= cLarger1 * (1-expm1a) * (expPhiApow2k / expm1a);
        }

        if (cSmaller2 > 0) {
            f  -= cSmaller2 * std::log(expm1bx);
            double tmp = cSmaller2 * (1-expm1bx) * (expPhiBXpow2k / expm1bx);
            fb -= tmp * expPhiBdivBX;
            fx -= tmp * expPhiXdivBX;
        }

        if (cLarger2 > 0) {
            f  -= cLarger2 * std::log(expm1b);
            fb -= cLarger2 * (1-expm1b) * (expPhiBpow2k / expm1b);
        }

        if (cEqual > 0) {

            double expm1x = -std::expm1(-expPhiXpow2k);
            double x = expm1a * expm1b * (1 - expm1x) + expm1x;
            f  -= cEqual * std::log(x);
            fa -= cEqual * ((((1 - expm1ax) * expm1b)/x) * expPhiApow2k);
            fb -= cEqual * ((((1 - expm1bx) * expm1a)/x) * expPhiBpow2k);
            fx -= cEqual * ((((1 - expm1x) *(1 - expm1a * expm1b))/x) * expPhiXpow2k);
        }

        if (k <= q) {
            expPhiApow2k += expPhiApow2k;
            expPhiBpow2k += expPhiBpow2k;
            expPhiXpow2k += expPhiXpow2k;
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


    const MaxLikelihoodEstimator estimator(jointStatistic.getP(), jointStatistic.getQ());

    const double cardinalityAX = estimator(jointStatistic.get1Counts());
    const double cardinalityBX = estimator(jointStatistic.get2Counts());
    const double cardinalityABX = estimator(jointStatistic.getMaxCounts());

    const double cardinalityMin = estimator(jointStatistic.getMinCounts());

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

    const gsl_multimin_fdfminimizer_type *T =  gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, 3); // 3 dimensions
    gsl_multimin_fdfminimizer_set(s, &my_func, phi, std::log(initalStepFactor), 0.1);

    numIterations = 0;
    do
    {
        numIterations++;
        int status = gsl_multimin_fdfminimizer_iterate(s);

        if (status) {
            if (status == GSL_ENOPROG) {
                iterationAborted = true;
                break;
            }
            std::cout << "error!" << std::endl;
            exit(-1);
        }

        const gsl_vector* currentPhi = gsl_multimin_fdfminimizer_x(s);
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

    cardinalityA = std::exp(gsl_vector_get(s->x, 0)) * m;
    cardinalityB = std::exp(gsl_vector_get(s->x, 1)) * m;
    cardinalityX = std::exp(gsl_vector_get(s->x, 2)) * m;

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(phi);

}

#endif // _CARDINALITY_ESTIMATION_HPP_
