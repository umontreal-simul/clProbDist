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
#ifndef CLPROBDIST_EXPONENTIALDIST_H
#define CLPROBDIST_EXPONENTIALDIST_H

/*! @file exponential.h
 *  @brief API of the exponential distribution
 *
 *  Implementation of \ref clProbDist_template.h for the exponential distribution,
 *  adapted from \cite iLEC08j .
 */

#include "clProbDist/clProbDist.h"
#include "clProbDist/continuous.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! @brief Exponential distribution object [**device**]
 *
 *  A structure that represents an exponential distribution object.
 */
typedef struct _clprobdistExponential clprobdistExponential;


/*! @name Functions to create and destroy distribution objects
 *
 * @{
 */

/*! @brief @copybrief clprobdistCreate()
 *
 *  Create a new exponential distribution object.
 *  Since this function allocates memory for the new distribution object;
 *  clprobdistDestroy() must be called to release the allocated memory.
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     New distribution object.
 */

clprobdistExponential* clprobdistExponentialCreate(cl_double lambda, size_t* bufSize, clprobdistStatus* err);

/*! @copydoc clprobdistDestroy()
 */
clprobdistStatus clprobdistExponentialDestroy(clprobdistExponential* dist);

/*! @} */


/*! @name Functions for use with a distribution object
 *
 * @{
 */

/*! @copydoc clprobdistDensityWithObject()
 */
cl_double clprobdistExponentialDensityWithObject(const clprobdistExponential* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistCDFWithObject()
 */
cl_double clprobdistExponentialCDFWithObject(const clprobdistExponential* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDFWithObject()
 */
cl_double clprobdistExponentialComplCDFWithObject(const clprobdistExponential* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDFWithObject()
 */
cl_double clprobdistExponentialInverseCDFWithObject(const clprobdistExponential* dist, cl_double u, clprobdistStatus* err);


/*! @copydoc clprobdistMeanWithObject()
 */
cl_double clprobdistExponentialMeanWithObject(const clprobdistExponential* dist, clprobdistStatus* err);

/*! @copydoc clprobdistVarianceWithObject()
 */
cl_double clprobdistExponentialVarianceWithObject(const clprobdistExponential* dist, clprobdistStatus* err);

/*! @copydoc clprobdistStdDeviationWithObject()
 */
cl_double clprobdistExponentialStdDeviationWithObject(const clprobdistExponential* dist, clprobdistStatus* err);

/*! @brief Return the value of the inverse mean \f$\lambda\f$  [**device**]
 */
cl_double clprobdistExponentialGetLambda(const clprobdistExponential* dist, clprobdistStatus* err);

/*! @brief Change the value of the inverse mean \f$\lambda\f$  [**device**]
 */
clprobdistStatus clprobdistExponentialSetLambda(clprobdistExponential* dist, cl_double newlambda);

/*! @} */


/*! @name Functions for use with explicit distribution parameters
 *
 *  @{
 */

/*! @brief @copybrief clprobdistDensity()
 *
 *  @see clprobdistExponentialDensityWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistExponentialDensity(cl_double lambda, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDF()
 *
 *  @see clprobdistExponentialCDFWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistExponentialCDF(cl_double lambda, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDF()
 *
 *  @see clprobdistExponentialComplCDFWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistExponentialComplCDF(cl_double lambda, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDF()
 *
 *  @see clprobdistExponentialInverseCDFWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
cl_double clprobdistExponentialInverseCDF(cl_double lambda, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMean()
 *
 *  @see clprobdistExponentialMeanWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Mean of the distribution.
 */
cl_double clprobdistExponentialMean(cl_double lambda, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVariance()
 *
 *  @see clprobdistExponentialVarianceWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistExponentialVariance(cl_double lambda, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviation()
 *
 *  @see clprobdistExponentialStdDeviationWithObject()
 *
 *  @param[in]  lambda  Value of the inverse mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistExponentialStdDeviation(cl_double lambda, clprobdistStatus* err);

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
