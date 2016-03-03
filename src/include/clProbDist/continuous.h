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
#ifndef CLPROBDIST_CONTINUOUSDIST_H
#define CLPROBDIST_CONTINUOUSDIST_H

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

/*! @cond PRIVATE */

/*! @file continuous.h
 *  @brief Common definitions for continuous distributions
 *
 *  Implementation of continuous distributions should include this header
 *
 *  This file provides default implementations for \f$\bar F(x)\f$
 *  and for \f$F^{-1}(u)\f$, the latter using the Brent-Dekker method (see
 *  \cite iBRE71a and \cite iBRE73a) to find the inverse of a generic distribution
 *  function \f$F(f)\f$.
 *  This is a robust root-finding method, appropriate to find a root \f$x\f$ of
 *  \f$F(x)=u\f$ when \f$u\f$ is known and if we have already a method to compute
 *  \f$F\f$.
 *  This is a special case.
 *
 */
typedef struct _clprobdistContinuousDist{
	int decPrec;             //Decimal degits of precision
	// [supportA, supportB] is the support of the pdf(x)
	cl_double supportA; //NEGATIVE INFINITY
	cl_double supportB; //POSITIVE INFINITY
} clprobdistContinuous;

// None of these functions are currently used by the library, but they could be
// used in the future when more distributions are implemented.

#if 0
/***********************************
* Constructor & destructor
***********************************/
cl_int clprobdistContinuousGetDecPrec(clprobdistContinuous* distObj, clprobdistStatus* err);

cl_double clprobdistContinuousGetXinf(clprobdistContinuous* distObj, clprobdistStatus* err);

cl_double clprobdistContinuousGetXsup(clprobdistContinuous* distObj, clprobdistStatus* err);

clprobdistStatus clprobdistContinuousSetXinf(cl_double xa, clprobdistContinuous* distObj);

clprobdistStatus clprobdistContinuousSetXsup(cl_double xb, clprobdistContinuous* distObj);



/***********************************
* Abstract functions
***********************************/
cl_double clprobdistContinuousCDF(cl_double x, clprobdistStatus* err);

clprobdistStatus findInterval(cl_double u, cl_double* iv, clprobdistContinuous* distObj, clprobdistStatus* err);

cl_double clprobdistContinuousDensity(cl_double x, clprobdistStatus* err);

cl_double clprobdistContinuousGetMean(clprobdistStatus* err);

cl_double clprobdistContinuousGetVariance(clprobdistStatus* err);

cl_double clprobdistContinuousGetStdDeviation(clprobdistStatus* err);

/***********************************
* Default implementations of functions
***********************************/
cl_double clprobdistContinuousComplCDF(cl_double x, clprobdistStatus* err);

cl_double clprobdistContinuousInverseBrent(cl_double a, cl_double b, cl_double u, cl_double tol, clprobdistContinuous* distObj, clprobdistStatus* err);

cl_double clprobdistContinuousInverseBisection(cl_double u, clprobdistContinuous* distObj, clprobdistStatus* err);

cl_double clprobdistContinuousInverseCDF(cl_double u, clprobdistContinuous* distObj, clprobdistStatus* err);
#endif

#endif /* CONTINUOUSDIST_H */
/*! @endcond */
