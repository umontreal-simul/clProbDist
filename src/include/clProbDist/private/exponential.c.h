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
#ifndef CLPROBDIST_PRIVATE_EXPONENTIALDIST_CH
#define CLPROBDIST_PRIVATE_EXPONENTIALDIST_CH

#ifndef _CLPROBDIST_EXPONENTIAL_OBJ_MEM
#define _CLPROBDIST_EXPONENTIAL_OBJ_MEM
#endif

struct _clprobdistExponential {
	clprobdistContinuous continuousDistObj;
	cl_double lambda;
};


constant const cl_double clprobdistExponentialXBIG = 100.0;
constant const cl_double clprobdistExponentialXBIGM = 1000.0;

clprobdistStatus clprobdistExponentialSetLambda(_CLPROBDIST_EXPONENTIAL_OBJ_MEM clprobdistExponential* distObj, cl_double newlambda)
{
	clprobdistStatus err = CLPROBDIST_SUCCESS;

	if (distObj->lambda <= 0)
		return clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);

	distObj->lambda = newlambda;
	distObj->continuousDistObj.supportA = 0.0;

	return err;
}

cl_double clprobdistExponentialDensity(cl_double lambda, cl_double x, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;

	if (lambda <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return x < 0 ? 0 : lambda * exp(-lambda * x);
}

cl_double clprobdistExponentialCDF(cl_double lambda, cl_double x, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;

	if (lambda <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	if (x <= 0.0)
		return 0.0;

	cl_double y = lambda * x;
	if (y >= clprobdistExponentialXBIG)
		return 1.0;

	return -expm1(-y);
}

cl_double clprobdistExponentialComplCDF(cl_double lambda, cl_double x, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;

	if (lambda <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}
	if (x <= 0.0)
		return 1.0;

	if (lambda*x >= clprobdistExponentialXBIGM)
		return 0.0;

	return exp(-lambda*x);
}

cl_double clprobdistExponentialInverseCDF(cl_double lambda, cl_double u, clprobdistStatus* err){

	if (err) *err = CLPROBDIST_SUCCESS;
	if (lambda <= 0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	if (u < 0.0 || u > 1.0)	{
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): u not in [0,1]", __func__);
		return -1;
	}

	if (u >= 1.0)
		return DBL_MAX;

	if (u <= 0.0)
		return 0.0;

	return -log1p(-u) / lambda; // log1p is defined on OpenCL C
}

cl_double clprobdistExponentialMean(cl_double lambda, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;
	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return (1 / lambda);
}

cl_double clprobdistExponentialVariance(cl_double lambda, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;
	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return (1 / (lambda * lambda));
}

cl_double clprobdistExponentialStdDeviation(cl_double lambda, clprobdistStatus* err)
{
	*err = CLPROBDIST_SUCCESS;
	if (lambda <= 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);
		return -1;
	}

	return (1 / lambda);
}


cl_double clprobdistExponentialGetLambda(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, clprobdistStatus* err)
{
	if (err) *err = CLPROBDIST_SUCCESS;
	return distObj->lambda;
}


cl_double clprobdistExponentialDensityWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, cl_double x, clprobdistStatus* err) {
	return clprobdistExponentialDensity(distObj->lambda, x, err);
}

cl_double clprobdistExponentialCDFWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, cl_double x, clprobdistStatus* err) {
	return clprobdistExponentialCDF(distObj->lambda, x, err);
}

cl_double clprobdistExponentialComplCDFWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, cl_double x, clprobdistStatus* err) {
	return clprobdistExponentialComplCDF(distObj->lambda, x, err);
}

cl_double clprobdistExponentialInverseCDFWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, cl_double u, clprobdistStatus* err) {
	return clprobdistExponentialInverseCDF(distObj->lambda, u, err);
}


cl_double clprobdistExponentialMeanWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, clprobdistStatus* err){
	  return clprobdistExponentialMean(distObj->lambda, err);
}

cl_double clprobdistExponentialVarianceWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, clprobdistStatus* err){
	  return clprobdistExponentialVariance(distObj->lambda, err);
}

cl_double clprobdistExponentialStdDeviationWithObject(_CLPROBDIST_EXPONENTIAL_OBJ_MEM const clprobdistExponential* distObj, clprobdistStatus* err){
	  return clprobdistExponentialStdDeviation(distObj->lambda, err);
}


#endif
