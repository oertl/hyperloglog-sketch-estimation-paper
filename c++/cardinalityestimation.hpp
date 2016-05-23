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


int getPFromNumberRegisters(int m) {
    int p;
    std::frexp(m, &p);
    p -= 1;
    assert(m == size_t(1) << p);
    return p;
}

int getPFromCounts(const std::vector<int>& c) {
    int m = std::accumulate(c.begin(), c.end(), 0);
    return getPFromNumberRegisters(m);
}

int getPFromJointStatistic(const std::vector<int>& jointStatistic) {
    size_t l = jointStatistic.size()/5;
    const int m = std::accumulate(&jointStatistic[l*0], &jointStatistic[l*1], 0) + std::accumulate(&jointStatistic[l*2], &jointStatistic[l*3], 0) + std::accumulate(&jointStatistic[l*3], &jointStatistic[l*4], 0);
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
            int e;
            frexp(x, &e);
            double xPrime = ldexp(x, -std::max(kMaxPrime+1, e+2));
            double xPrime2 = xPrime*xPrime;
            double h = xPrime - xPrime2/3 + (xPrime2*xPrime2)*(1./45. - xPrime2/472.5);
            for (int k = e; k >= kMaxPrime; --k) {
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
    return 1./(2.*std::log(2));
    /*assert(m >= 16);
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
    return alpha;*/
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

    static double sigma(int c, int m, int& numIterations) {
        numIterations = 0;
        if (c==m) return std::numeric_limits<double>::infinity();
        double powerOfX = static_cast<double>(c)/static_cast<double>(m);
        double result = powerOfX;
        double previousResult;
        double powerOf2 = 1;
        do {
            numIterations += 1;
            powerOfX *= powerOfX;
            previousResult = result;
            result += powerOfX * powerOf2;
            powerOf2 += powerOf2;
        } while(previousResult < result);
        return m*result;
    }

    static double tau(int c, int m, int& numIterations) {
        numIterations = 0;
        double result = static_cast<double>(c)/static_cast<double>(m);
        double powerOfX  = 1 - result;
        double previousResult;
        double powerOf2 = 0.5;
        do {
            numIterations += 1;
            powerOfX = std::sqrt(powerOfX);
            previousResult = result;
            result -= (1 - powerOfX)*powerOf2;
            powerOf2 *= 0.5;
        } while(previousResult > result);
        return m*result;
    }

    static std::vector<double> initSigma(int m) {
        int numIterations;
        std::vector<double> result(m+2);
        for (int c = 0; c < m+2; ++c) {
            result.push_back(sigma(c, m, numIterations));
        }
        return result;
    }

    static std::vector<double> initTau(int m) {
        int numIterations;
        std::vector<double> result(m+2);
        for (int c = 0; c < m+2; ++c) {
            result.push_back(tau(c, m, numIterations));
        }
        return result;
    }


public:

    CorrectedRawEstimator(const int p_, const int q_) :p(p_), q(q_) {}

    double estimate(const std::vector<int>& c, int& numSmallCorrectionIterations, int& numLargeCorrectionIterations) const {

        numSmallCorrectionIterations = 0;
        numLargeCorrectionIterations = 0;

        double z = 0.5 * tau(c[q+1], m, numLargeCorrectionIterations);
        for (int k = q; k >= 1; --k) {
            z += 0.5*c[k];
        }
        z += sigma(c[0], m, numSmallCorrectionIterations);
        return m2alpha/z;
    }

    double estimate2(const std::vector<int>& c) const {

        double z = 0.5 * tauValues[c[q+1]];
        for (int k = q; k >= 1; --k) {
            z += 0.5*c[k];
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
        if (e <= ldexp(1., 32)) { // is not part of original estimate
            eStar = -std::pow(2., 32)*std::log1p(-ldexp(e, -32));
        }
        else {
            eStar = std::numeric_limits<double>::infinity();
        }
    }
    return eStar;
}

// lower bound estimate obtained by applying Jensen's inequality
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

// using log(1+x) >= 2*x/(x+2)
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

void inclusionExclusionTwoHyperLogLogEstimation(const std::vector<int>& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX) {
    const int q = jointStatistic.size()/5-2;
    const int p = getPFromJointStatistic(jointStatistic);
    std::vector<int> countsAX(q+2);
    std::vector<int> countsBX(q+2);
    std::vector<int> countsABX(q+2);
    for (std::size_t i = 0; i < q+2; ++i) {
        countsAX[i]  = jointStatistic[(q+2)*0 + i] + jointStatistic[(q+2)*2 + i] + jointStatistic[(q+2)*4 + i];
        countsBX[i]  = jointStatistic[(q+2)*0 + i] + jointStatistic[(q+2)*1 + i] + jointStatistic[(q+2)*3 + i];
        countsABX[i] = jointStatistic[(q+2)*0 + i] + jointStatistic[(q+2)*1 + i] + jointStatistic[(q+2)*4 + i];
    }

    const MaxLikelihoodEstimator estimator(p,q);

    const double cardinalityAX = estimator(countsAX);
    const double cardinalityBX = estimator(countsBX);
    const double cardinalityABX = estimator(countsABX);

    cardinalityA = cardinalityABX - cardinalityBX;
    cardinalityB = cardinalityABX - cardinalityAX;
    cardinalityX = std::max(0., cardinalityBX + cardinalityAX - cardinalityABX);

}



// the minimum of this function gives the max likelihood estimate
void eval_joint_log_likelihood_function_and_derivatives(
    const std::vector<int>& jointStatistic,
    const double la,
    const double lb,
    const double lx,
    double& f,
    double& fa,
    double& fb,
    double& fx) {

    const int q = jointStatistic.size()/5-2;

    assert((q+2)*5==jointStatistic.size());

    f = 0.;
    fa = 0.;
    fb = 0.;
    fx = 0.;

    double term1a = 0;
    double term1b = 0;
    double term1x = 0;

    for (int k = q; k >= 0; --k) {
        double pow2k = std::ldexp(1., -k);
        term1a += (jointStatistic[(q+2)*0+k]+jointStatistic[(q+2)*2+k]+jointStatistic[(q+2)*4+k]) * pow2k;
        term1b += (jointStatistic[(q+2)*0+k]+jointStatistic[(q+2)*1+k]+jointStatistic[(q+2)*3+k]) * pow2k;
        term1x += (jointStatistic[(q+2)*0+k]+jointStatistic[(q+2)*2+k]+jointStatistic[(q+2)*3+k]) * pow2k;
    }

    for (int k = q+1; k >= 1; --k) {

        const int cCenter = jointStatistic[(q+2)*0+k];
        const int cUp     = jointStatistic[(q+2)*1+k];
        const int cRight  = jointStatistic[(q+2)*2+k];
        const int cDown   = jointStatistic[(q+2)*3+k];
        const int cLeft   = jointStatistic[(q+2)*4+k];

        double pow2k = std::ldexp(1., -std::min(q, k));

        double expm1a;
        double expm1b;
        double expm1ax;
        double expm1bx;

        if (cLeft > 0 || cCenter > 0) {
            expm1a = -std::expm1(-la * pow2k);
        }
        if (cUp > 0 || cCenter > 0) {
            expm1b = -std::expm1(-lb * pow2k);
        }
        if (cRight > 0 || cCenter > 0) {
            expm1ax = -std::expm1(-(la+lx) * pow2k);
        }
        if (cDown > 0 || cCenter > 0) {
            expm1bx = -std::expm1(-(lb+lx) * pow2k);
        }

        if (cRight > 0) {
            f  -= cRight * std::log(expm1ax);
            double tmp = cRight * pow2k * (1-expm1ax) / expm1ax;
            fa -= tmp;
            fx -= tmp;
        }

        if (cLeft > 0) {
            f  -= cLeft * std::log(expm1a);
            fa -= cLeft * pow2k * (1-expm1a) / expm1a;
        }

        if (cDown > 0) {
            f  -= cDown * std::log(expm1bx);
            double tmp = cDown * pow2k * (1-expm1bx) / expm1bx;
            fb -= tmp;
            fx -= tmp;
        }

        if (cUp > 0) {
            f  -= cUp * std::log(expm1b);
            fb -= cUp * pow2k * (1-expm1b) / expm1b;
        }

        if (cCenter > 0) {

            double expm1x = -std::expm1(-lx * pow2k);
            double x = expm1a * expm1b * (1 - expm1x) + expm1x;
            f  -= cCenter * std::log(x);
            fa -= cCenter * pow2k * ((1 - expm1ax) * expm1b)/x;
            fb -= cCenter * pow2k * ((1 - expm1bx) * expm1a)/x;
            fx -= cCenter * pow2k * ((1 - expm1x) *(1 - expm1a * expm1b))/x;
        }
    }

    f += term1a * la + term1b * lb + term1x * lx;
    fa += term1a;
    fb += term1b;
    fx += term1x;
}

void log_likelihood_function_value_and_derivatives_for_gsl(const gsl_vector *x, void *params, double *g, gsl_vector *dg)
{
    const double la = std::exp(gsl_vector_get(x, 0));
    const double lb = std::exp(gsl_vector_get(x, 1));
    const double lx = std::exp(gsl_vector_get(x, 2));
    const std::vector<int>& jointStatistic= *((std::vector<int>*)params);

    double f, fa, fb, fx;

    eval_joint_log_likelihood_function_and_derivatives(
        *((std::vector<int>*)params),
        la, lb, lx,
        f, fa, fb, fx);

    //if (f > std::numeric_limits<double>::max()) {
    //  f = std::numeric_limits<double>::max();
    //}
    *g = f;
    gsl_vector_set(dg, 0, fa*la);
    gsl_vector_set(dg, 1, fb*lb);
    gsl_vector_set(dg, 2, fx*lx);

}

void log_likelihood_function_derivatives_for_gsl(const gsl_vector *x, void *params, gsl_vector *df)
{
    double tmp;
    log_likelihood_function_value_and_derivatives_for_gsl(x,params, &tmp, df);
}

double log_likelihood_function_value_for_gsl(const gsl_vector *x, void *params)
{
    gsl_vector *grad = gsl_vector_alloc(3);
    double f;
    log_likelihood_function_value_and_derivatives_for_gsl(x, params, &f, grad);
    gsl_vector_free (grad);
    return f;
}


void maxLikelihoodTwoHyperLogLogEstimation(const std::vector<int>& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, bool& maxNumIterationsReached, bool& iterationAborted) {

    const int maxNumIterations = 10000;

    size_t l = jointStatistic.size()/5;
    const int m = std::accumulate(&jointStatistic[l*0], &jointStatistic[l*1], 0) + std::accumulate(&jointStatistic[l*2], &jointStatistic[l*3], 0) + std::accumulate(&jointStatistic[l*3], &jointStatistic[l*4], 0);

    double initialCardinalityA;
    double initialCardinalityB;
    double initialCardinalityX;
    inclusionExclusionTwoHyperLogLogEstimation(jointStatistic, initialCardinalityA, initialCardinalityB, initialCardinalityX);

    // set initial vector
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector_set (x, 0, std::log(std::max(1.,initialCardinalityA)/m));
    gsl_vector_set (x, 1, std::log(std::max(1.,initialCardinalityB)/m));
    gsl_vector_set (x, 2, std::log(std::max(1.,initialCardinalityX)/m));

    // set initial simplex size
    gsl_vector *step_size = gsl_vector_alloc(3);
    gsl_vector_set_all(step_size, 1.);


    gsl_multimin_function my_func;

    my_func.n = 3;
    my_func.f = log_likelihood_function_value_for_gsl;
    my_func.params = &const_cast<std::vector<int>&>(jointStatistic);


    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    const gsl_multimin_fminimizer_type *T =  gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 3); // 3 dimensions
    gsl_multimin_fminimizer_set (s, &my_func, x, step_size);

    int status;
    size_t iter = 0;
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) {
            if (status == GSL_ENOPROG) {
                iterationAborted = true;
                break;
            }
            std::cout << "error!" << std::endl;
            exit(-1);
        }
        //double cardmin = gsl_vector_min(gsl_multimin_fminimizer_x(s));
        double size = gsl_multimin_fminimizer_size(s);
        //std::cout << "size = " << size << std::endl;
        //std::cout << "cardmin = " << cardmin << std::endl;
        status = gsl_multimin_test_size(size, 1e-5);
        if (status != GSL_CONTINUE) {
            break;
        }
        if(iter >= maxNumIterations) {
            maxNumIterationsReached = true;
            break;
        }

    }
    while (true);


    //std::cout << iter << std::endl;

    cardinalityA = std::exp(gsl_vector_get(s->x, 0)) * m;
    cardinalityB = std::exp(gsl_vector_get(s->x, 1)) * m;
    cardinalityX = std::exp(gsl_vector_get(s->x, 2)) * m;

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(step_size);

}

void maxLikelihoodTwoHyperLogLogEstimation2(const std::vector<int>& jointStatistic, double& cardinalityA, double& cardinalityB, double& cardinalityX, bool& maxNumIterationsReached, bool& iterationAborted) {

    const int maxNumIterations = 10000;

    maxNumIterationsReached = false;
    iterationAborted = false;

    size_t l = jointStatistic.size()/5;
    const int m = std::accumulate(&jointStatistic[l*0], &jointStatistic[l*1], 0) + std::accumulate(&jointStatistic[l*2], &jointStatistic[l*3], 0) + std::accumulate(&jointStatistic[l*3], &jointStatistic[l*4], 0);

    double initialCardinalityA;
    double initialCardinalityB;
    double initialCardinalityX;
    inclusionExclusionTwoHyperLogLogEstimation(jointStatistic, initialCardinalityA, initialCardinalityB, initialCardinalityX);

    // set initial vector
    double lastA = std::log(std::max(1.,initialCardinalityA)/m);
    double lastB = std::log(std::max(1.,initialCardinalityB)/m);
    double lastX = std::log(std::max(1.,initialCardinalityX)/m);
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector_set (x, 0, lastA);
    gsl_vector_set (x, 1, lastB);
    gsl_vector_set (x, 2, lastX);

    gsl_multimin_function_fdf my_func;

    my_func.n = 3;
    my_func.f = log_likelihood_function_value_for_gsl;
    my_func.df = log_likelihood_function_derivatives_for_gsl;
    my_func.fdf = log_likelihood_function_value_and_derivatives_for_gsl;
    my_func.params = &const_cast<std::vector<int>&>(jointStatistic);

    const gsl_multimin_fdfminimizer_type *T =  gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 3); // 3 dimensions
    gsl_multimin_fdfminimizer_set (s, &my_func, x, 1, 0.1);

    size_t iter = 0;
    do
    {
        iter++;
        int status = gsl_multimin_fdfminimizer_iterate(s);

        if (status) {
            if (status == GSL_ENOPROG) {
                iterationAborted = true;
                break;
            }
            std::cout << "error!" << std::endl;
            exit(-1);
        }

        const gsl_vector* current_x = gsl_multimin_fdfminimizer_x(s);
        double currentA = gsl_vector_get(current_x, 0);
        double currentB = gsl_vector_get(current_x, 1);
        double currentX = gsl_vector_get(current_x, 2);

        if (
            std::fabs(currentA - lastA) <= 1e-5 &&
            std::fabs(currentB - lastB) <= 1e-5 &&
            std::fabs(currentX - lastX) <= 1e-5) {
            break;
        }

        lastA = currentA;
        lastB = currentB;
        lastX = currentX;

        if(iter >= maxNumIterations) {
            maxNumIterationsReached = true;
            break;
        }
    }
    while (true);

    cardinalityA = std::exp(gsl_vector_get(s->x, 0)) * m;
    cardinalityB = std::exp(gsl_vector_get(s->x, 1)) * m;
    cardinalityX = std::exp(gsl_vector_get(s->x, 2)) * m;

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

}



#endif // _CARDINALITY_ESTIMATION_HPP_
