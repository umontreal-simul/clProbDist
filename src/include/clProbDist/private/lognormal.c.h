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

#pragma once
#ifndef CLPROBDIST_PRIVATE_LOGNORMALDIST_CH
#define CLPROBDIST_PRIVATE_LOGNORMALDIST_CH

#ifndef _CLPROBDIST_LOGNORMAL_OBJ_MEM
#define _CLPROBDIST_LOGNORMAL_OBJ_MEM
#endif

#ifdef __cplusplus
extern "C"
{
#endif

struct _clprobdistlLognormal {
	clprobdistContinuous contDistObj;
	cl_double mu;
	cl_double sigma;
};

constant cl_double clprobdistLognormalLN2 = 0.6931471805599453;
constant cl_double clprobdistLognormalPI = 3.14159265358979323846;

constant const cl_double clprobdistLognormalXBIG = 100.0;
constant const cl_double clprobdistLognormalXBIGM = 1000.0;

//***********************************************************************
// Implementation of static functions
//***********************************************************************

cl_double clprobdistLognormalDensity(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err) {
	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	if (x <= 0)
		return 0;

	cl_double diff = log(x) - mu;

	return exp(-diff * diff / (2 * sigma * sigma)) / (sqrt(2 * clprobdistLognormalPI) * sigma * x);
}

cl_double clprobdistLognormalCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err) {
	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}

	if (x <= 0.0)
		return 0.0;
	return clprobdistStdNormalCDF((log(x) - mu) / sigma, err);
}
cl_double clprobdistLognormalComplCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err) {
	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}
	if (x <= 0.0)
		return 1.0;
	return clprobdistStdNormalComplCDF((log(x) - mu) / sigma, err);
}
cl_double clprobdistLognormalInverseCDF(cl_double mu, cl_double sigma, cl_double u, clprobdistStatus* err) {
	cl_double t, v;

	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}

	if (u > 1.0 || u < 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u is not in [0, 1]", __func__);
		return -1;
	}

	if (DBL_EPSILON >= 1.0 - u)
		return DBL_MAX; // Double.POSITIVE_INFINITY;

	if (u <= 0.0)
		return 0.0;

	t = clprobdistStdNormalInverseCDF(u, err); // NormalDist.inverseF01(u);
	v = mu + sigma * t;

	if ((t >= clprobdistLognormalXBIG) || (v >= DBL_MAX_EXP * clprobdistLognormalLN2))
		return DBL_MAX;
	if ((t <= -clprobdistLognormalXBIG) || (v <= -DBL_MAX_EXP * clprobdistLognormalLN2))
		return 0.0;

	return exp(v);
}
cl_double clprobdistLognormalMean(cl_double mu, cl_double sigma, clprobdistStatus* err) {
	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}

	return (exp(mu + (sigma * sigma) / 2.0));
}
cl_double clprobdistLognormalVariance(cl_double mu, cl_double sigma, clprobdistStatus* err) {
	if (sigma <= 0.0) {
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): sigma <= 0", __func__);
		return 0;
	}

	return (exp(2.0 * mu + sigma * sigma) * (exp(sigma * sigma) - 1.0));
}
cl_double clprobdistLognormalStdDeviation(cl_double mu, cl_double sigma, clprobdistStatus* err) {

	return sqrt(clprobdistLognormalVariance(mu, sigma, err));
}

//***********************************************************************
// Implementation with Dist object
//***********************************************************************

cl_double clprobdistLognormalDensityWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err) {
	return clprobdistLognormalDensity(dist->mu, dist->sigma, x, err);
}

cl_double clprobdistLognormalCDFWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err) {
	return clprobdistLognormalCDF(dist->mu, dist->sigma, x, err);
}

cl_double clprobdistLognormalComplCDFWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err) {
	return clprobdistLognormalComplCDF(dist->mu, dist->sigma, x, err);
}

cl_double clprobdistLognormalInverseCDFWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, cl_double u, clprobdistStatus* err) {
	return clprobdistLognormalInverseCDF(dist->mu, dist->sigma, u, err);
}

cl_double clprobdistLognormalMeanWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, clprobdistStatus* err) {
	return clprobdistLognormalMean(dist->mu, dist->sigma, err);
}

cl_double clprobdistLognormalVarianceWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, clprobdistStatus* err) {
	return clprobdistLognormalVariance(dist->mu, dist->sigma, err);
}

cl_double clprobdistLognormalStdDeviationWithObject(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, clprobdistStatus* err) {
	return clprobdistLognormalStdDeviation(dist->mu, dist->sigma, err);
}

cl_double clprobdistLognormalGetMu(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return dist->mu;
}

cl_double clprobdistLognormalGetSigma(_CLPROBDIST_LOGNORMAL_OBJ_MEM const clprobdistLognormal* dist, clprobdistStatus* err) {
	if (err) *err = CLPROBDIST_SUCCESS;
	return dist->sigma;
}

#ifdef __cplusplus
}
#endif

#endif
