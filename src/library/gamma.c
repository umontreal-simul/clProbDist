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
#include <float.h>

#include "../include/clProbDist/clProbDist.h"
#include "private.h"

#include "../include/clProbDist/continuous.h"
#include "../include/clProbDist/exponential.h"
#include "../include/clProbDist/normal.h"
#include "../include/clProbDist/gamma.h"

// x infinity for some distributions
//extern const cl_double clprobdistXBIG;
//extern const cl_double clprobdistXBIGM;

// clprobdistEPSARRAY[j]: Epsilon required for j decimal degits of precision
//extern const cl_double clprobdistEPSARRAY[];

// code that is common to host and device
#include "../include/clProbDist/private/gamma.c.h"

clprobdistGamma* clprobdistGammaCreate(cl_double alpha, cl_double lambda, int decprec, size_t* bufSize, clprobdistStatus* err)
{
	clprobdistStatus err_ = CLPROBDIST_SUCCESS;

	//Check params
	if (alpha <= 0.0)
		err_ = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): alpha <= 0", __func__);

	if (decprec <= 0)
		err_ = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): d <= 0", __func__);

	if (lambda <= 0)
		err_ = clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "%s(): lambda <= 0", __func__);

	if (err_ != CLPROBDIST_SUCCESS) {
	    if (err) *err = err_;
		return NULL;
	}

	clprobdistGamma* gammaDist = (clprobdistGamma*)malloc(sizeof(clprobdistGamma));

	if (!gammaDist) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	if (bufSize)
		  *bufSize = sizeof(clprobdistGamma);
	if (err)
		  *err = CLPROBDIST_SUCCESS;

	gammaDist->contDistObj.decPrec = decprec;
	gammaDist->contDistObj.supportA = 0.0;
	gammaDist->contDistObj.supportB = DBL_MAX;

	gammaDist->params.alpha = alpha;
	gammaDist->params.lambda = lambda;
	gammaDist->params.logFactor = alpha * log(lambda) - lgamma(alpha);

	return gammaDist;
}
clprobdistStatus clprobdistGammaDestroy(clprobdistGamma* distObj){
	if (distObj != NULL)
		free(distObj);

	return CLPROBDIST_SUCCESS;
}



