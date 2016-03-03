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
#include "../include/clProbDist/exponential.h"

// x infinity for some distributions
//extern const cl_double clprobdistXBIG;
//extern const cl_double clprobdistXBIGM;

// code that is common to host and device
#include "../include/clProbDist/private/exponential.c.h"

clprobdistExponential* clprobdistExponentialCreate(cl_double lambda, size_t* bufSize, clprobdistStatus* err)
{
	clprobdistExponential* expoDistObj = (clprobdistExponential*)malloc(sizeof(clprobdistExponential));

	if (!expoDistObj) {
		  if (err)
			    *err = clprobdistSetErrorString(CLPROBDIST_OUT_OF_RESOURCES, "%s(): could not allocate memory for distribution", __func__);
		  return NULL;
	}

	if (bufSize)
		  *bufSize = sizeof(clprobdistExponential);
	if (err)
		  *err = CLPROBDIST_SUCCESS;

	expoDistObj->continuousDistObj.decPrec = 15;
	expoDistObj->continuousDistObj.supportA = 0.0;
	expoDistObj->continuousDistObj.supportB = DBL_MAX;
	expoDistObj->lambda = lambda;

	return expoDistObj;
}
clprobdistStatus clprobdistExponentialDestroy(clprobdistExponential* distObj)
{
	if (distObj != NULL)
		free(distObj);

	return CLPROBDIST_SUCCESS;
}





