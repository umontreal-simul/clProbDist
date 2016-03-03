/* This file is part of clProbDist.
 *
 * Copyright 2015-2016  Pierre L'Ecuyer, Universite de Montreal and Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Authors:
 *
 *   Nabil Kemerchou <kemerchn@iro.umontreal.ca>    (2015)
 *   David Munger <mungerd@iro.umontreal.ca>        (2015)
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 */

/*! @file clprobdistPoisson.c.h
*  @brief Code for the clprobdistPoisson generator common to the host and device
*/

#pragma once
#ifndef CLPROBDIST_PRIVATE_POISSONDIST_CH
#define CLPROBDIST_PRIVATE_POISSONDIST_CH

#ifndef _CLPROBDIST_POISSON_OBJ_MEM
#define _CLPROBDIST_POISSON_OBJ_MEM
#endif

#ifndef __OPENCL_C_VERSION__
#include <limits.h>
#endif

typedef struct _clprobdistPoissonParams {
	cl_int xmin;           // pdf[x] < EPSILON for x < xmin
	cl_int xmax;           // pdf[x] < EPSILON for x > xmax

	// xmed is such that cdf[xmed] >= 0.5 and cdf[xmed - 1] < 0.5.
	cl_int xmed;           // cdf[x] = F(x) for x <= xmed, and cdf[x] = bar_F(x) for x > xmed

	cl_int supportA;
	cl_int supportB;

	cl_double lambda;
	cl_int len;  // length of the CDF and PDF arrays

} clprobdistPoissonParams;

struct _clprobdistPoisson {
	clprobdistPoissonParams params;
	cl_double pdf[1];     // probability terms or mass distribution
	cl_double cdf[1];     // cumulative probabilities
};

constant const cl_double clprobdistPoisson_EPSILON = 1.0e-16;

cl_double clprobdistPoissonCDF_1(cl_double lambda, cl_int x, clprobdistStatus* err);

static cl_double clprobdistPoisson_CDF(cl_int x, _CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist){
	return dist->cdf[x + (dist->params.len - 1)];
}

static cl_double clprobdistPoisson_factorial(cl_int n, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (n < 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): n < 0", __func__);
		return -1;
	}
		
	cl_double T = 1;
	for (cl_int j = 2; j <= n; j++)
		T *= j;
	return T;
}

cl_double clprobdistPoissonProb(cl_double lambda, cl_int x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (x < 0) 
		return 0.0;

	if (lambda >= 100.0) {
		if ((cl_double)x >= 10.0*lambda)
			return 0.0;
	}
	else if (lambda >= 3.0) {
		if ((cl_double)x >= 100.0*lambda)
			return 0.0;			
	}
	else {
		if ((cl_double)x >= 200.0*fmax(1.0, lambda))
			return 0.0;
	}

	cl_double lambdaLIM = 20.0;
	cl_double Res;
	if (lambda < lambdaLIM && x <= 100)
		Res = exp(-lambda)*pow(lambda, x) / clprobdistPoisson_factorial(x, err);
	else {
		cl_double y = x*log(lambda) - lgamma(x + 1.0) - lambda;
		Res = exp(y);
	}
	return Res;
}

cl_double clprobdistPoissonCDF(cl_double lambda, cl_int x, clprobdistStatus* err);

cl_double clprobdistPoissonComplCDF(cl_double lambda, cl_int x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	//check Params
	if (lambda < 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	if (x <= 0) 
		return 1.0;

	if (lambda >= 100.0) {
		if ((cl_double)x >= 10.0*lambda)
			return 0.0;
	}
	else {
		if ((cl_double)x >= 100 + 100.0 * fmax(1.0, lambda))
			return 0.0;		
	}

	/* If lambda > lambdaLIM, we use the Chi2 distribution according to the
	exact relation, with 2x + 2 degrees of freedom

	cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*lambda)

	which also equals   1 - clprobdistGamma.cdf (x + 1, lambda) */

	cl_double lambdaLIM = 200.0;
	if (lambda > lambdaLIM)
		return clprobdistGammaCDF_1((cl_double)x, 15, lambda, err);


	if (x <= lambda)
		return 1.0 - clprobdistPoissonCDF_1(lambda, x - 1, err);
	

	// Naive computation: sum all prob. from i = x to i = oo
	cl_double term, sum;
	cl_int IMAX = 20;

	// Sum at least IMAX prob. terms from i = s to i = oo
	sum = term = clprobdistPoissonProb(lambda, x, err);
	
	cl_int i = x + 1;
	while (term > clprobdistPoisson_EPSILON || i <= x + IMAX) {
		term *= lambda / i;
		sum += term;
		i++;
	}

	return sum;
}

cl_double clprobdistPoissonCDF(cl_double lambda, cl_int x, clprobdistStatus* err) {
	cl_double ret = clprobdistPoissonCDF_1(lambda, x, err);
	if (ret < 0.0)
		return 1.0 - clprobdistPoissonComplCDF(lambda, x + 1, err);
	return ret;
}

cl_double clprobdistPoissonCDF_1(cl_double lambda, cl_int x, clprobdistStatus* err) {
	/*
	* On our machine, computing a value using gamma is faster than the
	* naive computation for dist->params.lambdalim > 200.0, slower for dist->params.lambdalim < 200.0
	*/
	if (err) *err = CLPROBDIST_SUCCESS;
	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	if (lambda == 0.0)
		return 1.0;
		
	if (x < 0)
		return 0.0;
		

	if (lambda >= 100.0) {
		if ((cl_double)x >= 10.0*lambda)
			return 1.0;		
	}
	else {
		if ((cl_double)x >= 100.0*fmax(1.0, lambda))
			return 1.0;
	}

	/* If lambda > lambdaLIM, use the Chi2 distribution according to the
	exact relation, with 2x + 2 degrees of freedom

	poisson (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

	which also equals 1 - gamma (x + 1, lambda) */

	cl_double lambdaLIM = 200.0;
	if (lambda > lambdaLIM)
		return clprobdistGammaComplCDF_1(x + 1.0, 15, lambda, err);

	if (x >= lambda)
		// this is a workaround to avoid birecusivity on GPU's
		// by returning -1.0, we let the caller know it must call ComplCDF()
		return -1.0;

	// Naive computation: sum all prob. from i = x
	cl_double sum = 1;
	cl_double term = 1;
	for (cl_int j = 1; j <= x; j++) {
		term *= lambda / j;
		sum += term;
	}
	return sum*exp(-lambda);
}

cl_int clprobdistPoissonInverseCDF(cl_double lambda, cl_double u, clprobdistStatus* err) {
	//Check Params
	if (u < 0.0 || u > 1.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u is not in range [0,1]", __func__);
		return -1;
	}

	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	if (u >= 1.0)
		return INT_MAX;

	if (u <= clprobdistPoissonProb(lambda, 0, err)) 
		return 0;

	//Calculate InverseCDF
	cl_int i;
	cl_double lambdaLIM = 700.0;

	if (lambda < lambdaLIM)
	{
		cl_double sumprev = -1.0;
		cl_double term = exp(-lambda);
		cl_double sum = term;
		i = 0;
		while (sum < u && sum > sumprev) {
			i++;
			term *= lambda / i;
			sumprev = sum;
			sum += term;
		}

		return i;

	}
	else {
		// When lambda is very large, the probabilities are empirically
		// negligible when i is far from lambda.  We start with a binary search
		// over [0,lambda] for the smallest value of i for which the probability
		// is non-negligible.
		// Then we perform the sequential search starting from that value of i.

		// Search for value of i smaller than u yet close enough to u.
		i = (cl_int)lambda;
		cl_double term = clprobdistPoissonProb(lambda, i, err);

		while ((term >= u) && (term > DBL_MIN )) {
			i /= 2;
			term = clprobdistPoissonProb(lambda, i, err);
		}

		if (term <= DBL_MIN ) {
			i *= 2;
			term = clprobdistPoissonProb(lambda, i, err);
			while (term >= u && (term > DBL_MIN )) {
				term *= i / lambda;
				i--;
			}
		}

		cl_int mid = i;
		cl_double sum = term;
		cl_double termid = term;

		// Begin the search here, at i = mid.
		// We start by computing F(i), i.e., by summing all non-negligible
		// probabilities for i <= mid.
		// Probabilities below clprobdistPoisson_EPSILON*u are deemed
		// negligible.
		while (term >= clprobdistPoisson_EPSILON*u && i > 0) {
			term *= i / lambda;
			sum += term;
			i--;
		}

		// Sequential search for smallest i such that F(i) >= u.
		term = termid;
		i = mid;
		cl_double prev = -1;
		if (sum < u) {
			while ((sum < u) && (sum > prev)) {
				i++;
				term *= lambda / i;
				prev = sum;
				sum += term;
			}
		}
		else {
			// The computed CDF is too big so we substract from it.
			sum -= term;
			while (sum >= u) {
				term *= i / lambda;
				i--;
				sum -= term;
			}
		}

		return i;
	}
}

cl_double clprobdistPoissonMean(cl_double lambda, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	//Check Params
	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	return lambda;
}

cl_double clprobdistPoissonVariance(cl_double lambda, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	//Check Params
	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	return lambda;
}

cl_double clprobdistPoissonStdDeviation(cl_double lambda, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	//Check Params
	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return -1;
	}

	return sqrt(lambda);
}


//Dynamic Functions

cl_double clprobdistPoissonProbWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	if (x < 0) 
		return 0.0;

	if (!dist->pdf)
		return clprobdistPoissonProb(dist->params.lambda, x, err);

	if (x > dist->params.xmax || x < dist->params.xmin)
		return clprobdistPoissonProb(dist->params.lambda, x, err);

	return dist->pdf[x - dist->params.xmin];
}

cl_double clprobdistPoissonCDFWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	cl_double Sum = 0.0;

	if (x < 0)
		return 0.0;
		
	if (dist->params.lambda == 0.0)
		return 1.0;

	/* For large lambda, we use the Chi2 distribution according to the exact
	relation, with 2x + 2 degrees of freedom

	cdf (lambda, x) = 1 - chiSquare (2x + 2, 2*lambda)

	which equals also 1 - gamma (x + 1, dist->params.lambda) */
	
	if (!dist->cdf)
		return clprobdistGammaComplCDF_1(x + 1.0, 15, dist->params.lambda, err);

	if (x >= dist->params.xmax)
		return 1.0;

	if (x < dist->params.xmin) {
		// Sum a few terms to get a few decimals far in the lower tail. One
		// could also call clprobdistGamma.barF instead.
		cl_int RMAX = 20;
		cl_int i;
		cl_double term; 
		Sum = term = clprobdistPoissonProb(dist->params.lambda, x, err);
		i = x;
		while (i > 0 && i >= x - RMAX) {
			term = term * i / dist->params.lambda;
			i--;
			Sum += term;
		}
		return Sum;
	}

	if (x <= dist->params.xmed)
		return clprobdistPoisson_CDF(x - dist->params.xmin, dist);
	else
		// We keep the complementary distribution in the upper part of cdf
		return 1.0 - clprobdistPoisson_CDF(x + 1 - dist->params.xmin, dist);
}

cl_double clprobdistPoissonComplCDFWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err) {
	/*
	* barF (lambda, x) = 1 - cdf (lambda, x - 1)
	*/
	if (err) *err = CLPROBDIST_SUCCESS;

	if (x <= 0) 
		return 1.0;

	/* For large dist->params.lambda,  we use the Chi2 distribution according to the exact
	relation, with 2x + 2 degrees of freedom

	cdf (lambda, x) = 1 - ChiSquare.cdf (2x + 2, 2*dist->params.lambda)
	cdf (lambda, x) = 1 - clprobdistGamma.cdf (x + 1, dist->params.lambda)
	*/

	if (!dist->cdf)
		return clprobdistGammaCDF_1((cl_double)x, 15, dist->params.lambda, err);

	if (x > dist->params.xmax)
		return clprobdistPoissonComplCDF(dist->params.lambda, x, err);

	if (x <= dist->params.xmin)
		return 1.0;

	if (x > dist->params.xmed){
		// We keep the complementary distribution in the upper part of cdf
		return clprobdistPoisson_CDF(x - dist->params.xmin, dist);
	}
	else
		return 1.0 - clprobdistPoisson_CDF(x - 1 - dist->params.xmin, dist);
}

cl_int clprobdistDiscreteInverseCDF(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, cl_double u, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;

	cl_int i=0, j=0, k=0;

	if (u < 0.0 || u > 1.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u is not in [0,1]", __func__);
		return -1;
	}
		
	if (u <= 0.0)
		return dist->params.supportA;
		
	if (u >= 1.0)
		return dist->params.supportB;
		

	// Remember: the upper part of cdf contains the complementary distribu-
	// tion for xmed < s <= xmax, and the lower part of cdf the
	// distribution for xmin <= x <= xmed

	if (u <= clprobdistPoisson_CDF(dist->params.xmed - dist->params.xmin, dist)) {
		// In the lower part of cdf
		if (u <= clprobdistPoisson_CDF(0, dist))
			return dist->params.xmin;
			
		i = 0;
		j = dist->params.xmed - dist->params.xmin;
		while (i < j) {
			k = (i + j) / 2;
			if (u > clprobdistPoisson_CDF(k, dist))
				i = k + 1;
			else
				j = k;
		}
	}
	else {
		// In the upper part of cdf
		u = 1 - u;
		if (u < clprobdistPoisson_CDF(dist->params.xmax - dist->params.xmin, dist))
			return dist->params.xmax;

			i = dist->params.xmed - dist->params.xmin + 1;
			j = dist->params.xmax - dist->params.xmin;
			while (i < j) {
				k = (i + j) / 2;
				if (u < clprobdistPoisson_CDF(k, dist))
					i = k + 1;
				else
					j = k;
			}
			i--;
		
	}
	return i + dist->params.xmin;
}

cl_int clprobdistPoissonInverseCDFWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, cl_double u, clprobdistStatus* err) {
	if ((!dist->cdf) || (u <= clprobdistPoisson_EPSILON))
		return clprobdistPoissonInverseCDF(dist->params.lambda, u, err);

	return clprobdistDiscreteInverseCDF(dist, u, err);

}

cl_double clprobdistPoissonMeanWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, clprobdistStatus* err) {
	return clprobdistPoissonMean(dist->params.lambda, err);
}

cl_double clprobdistPoissonVarianceWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, clprobdistStatus* err) {
	return clprobdistPoissonVariance(dist->params.lambda, err);
}

cl_double clprobdistPoissonStdDeviationWithObject(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, clprobdistStatus* err) {
	return clprobdistPoissonStdDeviation(dist->params.lambda, err);
}

cl_double clprobdistPoissonGetLambda(_CLPROBDIST_POISSON_OBJ_MEM const clprobdistPoisson* dist, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;
	return dist->params.lambda;
}


#endif // PRIVATE_POISSONDIST_CH
