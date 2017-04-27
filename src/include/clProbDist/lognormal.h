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
#ifndef CLPROBDIST_LOGNORMALDIST_H
#define CLPROBDIST_LOGNORMALDIST_H

/*! @file lognormal.h
*  @brief API of the lognormal distribution
*
*  Implementation of \ref clProbDist_template.h for the lognormal distribution,
*  adapted from \cite iLEC08j .
*/


#include "clProbDist/clProbDist.h"
#include "clProbDist/continuous.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! @brief Lognormal distribution object [**device**]
*
*  A structure that represents a lognormal distribution object.
*/
typedef struct _clprobdistlLognormal clprobdistLognormal;


//***********************************************************************
// Constructor and Destructor
//***********************************************************************

/*! @name Functions to create and destroy distribution objects
*
* @{
*/

/*! @brief @copybrief clprobdistCreate()
*
*  Create a new lognormal distribution object.
*  Since this function allocates memory for the new distribution object;
*  clprobdistDestroy() must be called to release the allocated memory.
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     New distribution object.
*/
clprobdistLognormal* clprobdistLognormalCreate(cl_double mu, cl_double sigma, size_t* bufSize, clprobdistStatus* err);

/*! @copydoc clprobdistDestroy()
*/
clprobdistStatus clprobdistLognormalDestroy(clprobdistLognormal* dist);

/*! @} */

//***********************************************************************
// Implementation with Dist object
//***********************************************************************

/*! @name Functions for use with a distribution object
*
* @{
*/

/*! @copydoc clprobdistDensityWithObject()
*/
cl_double clprobdistLognormalDensityWithObject(const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistCDFWithObject()
*/
cl_double clprobdistLognormalCDFWithObject(const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDFWithObject()
*/
cl_double clprobdistLognormalComplCDFWithObject(const clprobdistLognormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDFWithObject()
*/
cl_double clprobdistLognormalInverseCDFWithObject(const clprobdistLognormal* dist, cl_double u, clprobdistStatus* err);

/*! @copydoc clprobdistMeanWithObject()
*/
cl_double clprobdistLognormalMeanWithObject(const clprobdistLognormal* dist, clprobdistStatus* err);

/*! @copydoc clprobdistVarianceWithObject()
*/
cl_double clprobdistLognormalVarianceWithObject(const clprobdistLognormal* dist, clprobdistStatus* err);

/*! @copydoc clprobdistStdDeviationWithObject()
*/
cl_double clprobdistLognormalStdDeviationWithObject(const clprobdistLognormal* dist, clprobdistStatus* err);

/*! @brief Return the scale parameter.
*
*  @param[in]  dist    Distribution object.
*  @param[out] err     Error status variable, or \c NULL.
*  @return Scale parameter \f$\mu\f$.
*/
cl_double clprobdistLognormalGetMu(const clprobdistLognormal* dist, clprobdistStatus* err);

/*! @brief Return the shape parameter.
*
*  @param[in]  dist    Distribution object.
*  @param[out] err     Error status variable, or \c NULL.
*  @return Shape parameter \f$\sigma\f$.
*/
cl_double clprobdistLognormalGetSigma(const clprobdistLognormal* dist, clprobdistStatus* err);

/*! @} */


//***********************************************************************
// Implementation of static functions
//***********************************************************************

/*! @name Functions for use with explicit distribution parameters
*
*  @{
*/

/*! @brief @copybrief clprobdistDensity()
*
*  @see clprobdistLognormalDensityWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[in]  x       Value of \f$x\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Value of \f$f(x)\f$.
*/
cl_double clprobdistLognormalDensity(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDF()
*
*  @see clprobdistLognormalCDFWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[in]  x       Value of \f$x\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Value of \f$F(x)\f$.
*/
cl_double clprobdistLognormalCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDF()
*
*  @see clprobdistLognormalComplCDFWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[in]  x       Value of \f$x\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Value of \f$\bar F(x)\f$.
*/
cl_double clprobdistLognormalComplCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDF()
*
*  @see clprobdistLognormalInverseCDFWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[in]  u       Value of \f$u \in [0,1]\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Value of \f$F^{-1}(u)\f$.
*/
cl_double clprobdistLognormalInverseCDF(cl_double mu, cl_double sigma, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMean()
*
*  @see clprobdistLognormalMeanWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return	Mean of the distribution.
*/
cl_double clprobdistLognormalMean(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVariance()
*
*  @see clprobdistLognormalVarianceWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Variance of the distribution.
*/
cl_double clprobdistLognormalVariance(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviation()
*
*  @see clprobdistLognormalStdDeviationWithObject()
*
*  @param[in]  mu	Value of the scale parameter \f$\mu\f$.
*  @param[in]  sigma	Value of the shape parameter \f$\sigma\f$.
*  @param[out] err     Error status variable, or \c NULL.
*  @return     Standard deviation of the distribution.
*/
cl_double clprobdistLognormalStdDeviation(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
