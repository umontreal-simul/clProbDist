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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../include/clProbDist/clProbDist.h"
#include "private.h"

#include "../include/clProbDist/gamma.h"
#include "../include/clProbDist/poisson.h"

// code that is common to host and device
#include "../include/clProbDist/private/poisson.c.h"


const cl_double MAXlambda = 100000;
const cl_double EPS_EXTRA = 1.0e-6;

/***********************************
* Constructor & destructor
***********************************/

clprobdistPoisson* clprobdistPoissonInitialize(cl_double lambda, size_t *bufSize, clprobdistStatus* err) {

	//Allocate distribution parameters
	clprobdistPoissonParams* params = (clprobdistPoissonParams*)malloc(sizeof(clprobdistPoissonParams));

	if (!params) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	params->supportA = 0;
	params->supportB = INT_MAX;

	if (lambda < 0.0){
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda < 0", __func__);
		return NULL;
	}

	params->lambda = lambda;


	// For lambda > MAXlambda, we do not use pre-computed arrays
	if (lambda > MAXlambda) {
		//dist->pdf = NULL;
		//dist->cdf = NULL;
		if (err) *err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda > MAXlambda ", __func__);
		return NULL;
	}

	cl_double epsilon;
	int i, mid, Nmax;
	int imin, imax;
	cl_double sum;
	cl_double* P;    // Poisson probability terms
	cl_double* F;    // Poisson cumulative probabilities

	// In theory, the Poisson distribution has an infinite range. But
	// for i > Nmax, probabilities should be extremely small.
	Nmax = (int)(lambda + 16 * (2 + sqrt(lambda)));
	P = (cl_double*)malloc((1 + Nmax)*sizeof(cl_double));
	if (!P) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	mid = (int)lambda;
	epsilon = clprobdistPoisson_EPSILON * EPS_EXTRA / clprobdistPoissonProb(lambda, mid, err);

	// For large lambda, mass will lose a few digits of precision
	// We shall normalize by explicitly summing all terms >= epsilon
	sum = P[mid] = 1.0;

	// Start from the maximum and compute terms > epsilon on each side.
	i = mid;
	while (i > 0 && P[i] > epsilon) {
		P[i - 1] = P[i] * i / lambda;
		i--;
		sum += P[i];
	}
	params->xmin = imin = i;

	i = mid;
	while (P[i] > epsilon) {
		P[i + 1] = P[i] * lambda / (i + 1);
		i++;
		sum += P[i];
		if (i >= Nmax - 1) {
			Nmax *= 2;
			cl_double* nT = (cl_double*)malloc((1 + Nmax) * sizeof(cl_double));
			if (!nT) {
			  if (err)
				*err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
			  return NULL;
			}
			arraycopy(P, 0, nT, 0, (1 + Nmax) /*P.length*/);
			P = nT;
		}
	}
	params->xmax = imax = i;
	F = (cl_double*)malloc((1 + Nmax)*sizeof(cl_double));
	if (!F) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	// Renormalize the sum of probabilities to 1
	for (i = imin; i <= imax; i++)
		P[i] /= sum;

	// Compute the cumulative probabilities until F >= 0.5, and keep them in
	// the lower part of array, i.e. F[s] contains all P[i] for i <= s
	F[imin] = P[imin];
	i = imin;
	while (i < imax && F[i] < 0.5) {
		i++;
		F[i] = P[i] + F[i - 1];
	}
	// This is the boundary between F and 1 - F in the CDF
	params->xmed = i;

	// Compute the cumulative probabilities of the complementary distribution
	// and keep them in the upper part of the array. i.e. F[s] contains all
	// P[i] for i >= s
	F[imax] = P[imax];
	i = imax - 1;
	do {
		F[i] = P[i] + F[i + 1];
		i--;
	} while (i > params->xmed);

	/* Reset imin because we lose too much precision for a few terms near
	imin when we stop adding terms < epsilon. */
	i = imin;
	while (i < params->xmed && F[i] < clprobdistPoisson_EPSILON)
		i++;
	params->xmin = imin = i;

	/* Same thing with imax */
	i = imax;
	while (i > params->xmed && F[i] < clprobdistPoisson_EPSILON)
		i--;
	params->xmax = imax = i;

	// Create the Poisson  Distribution with flex array :
	params->len = (imax + 1 - imin);
	size_t bufSize_ = sizeof(clprobdistPoisson) + 2 * ((params->len - 1) * sizeof(cl_double));
	clprobdistPoisson* dist = (clprobdistPoisson*)malloc(bufSize_);
	if (!dist) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	//Set dist values:
	dist->params = *params;

	//Copy CDF and PDF values
	for (int i = 0; i < params->len; i++){
		dist->pdf[i] = P[imin + i];
		dist->cdf[i + (params->len - 1)] = F[imin + i];
	}

	/*int len = params->len;
	for (unsigned ix = 0; ix < len; ix++)
	printf("pdf %d , %20.19f \n", ix, dist->pdf[ix]);

	for (unsigned ix = 0; ix < len; ix++)
	printf("cdf %d , %20.19f \n", ix, dist->cdf[ix + (len - 1)]);*/

	if (bufSize)
		  *bufSize = bufSize_;
	if (err)
		  *err = CLPROBDIST_SUCCESS;

	return dist;
}

clprobdistPoisson* clprobdistPoissonCreate(cl_double lambda, size_t* bufSize, clprobdistStatus* err)
{
	//check params
	if (lambda <= 0.0){
	  if (err)
		*err = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);
	  return NULL;
	}

	//Initialize the params and histogrammes
	clprobdistPoisson* dist = clprobdistPoissonInitialize(lambda, bufSize, err);

	return dist;
}
clprobdistStatus clprobdistPoissonDestroy(clprobdistPoisson* dist)
{
	if (dist != NULL)
		free(dist);

	return CLPROBDIST_SUCCESS;
}

