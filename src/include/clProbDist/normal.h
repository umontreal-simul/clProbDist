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
#ifndef CLPROBDIST_NORMALDIST_H
#define CLPROBDIST_NORMALDIST_H

/*! @file normal.h
 *  @brief API of the normal distribution
 *
 *  Implementation of \ref clProbDist_template.h for the normal distribution,
 *  adapted from \cite iLEC08j .
 */

#include "clProbDist/clProbDist.h"
#include "clProbDist/continuous.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! @brief Normal distribution object [**device**]
 *
 *  A structure that represents a normal distribution object.
 */
typedef struct _clprobdistNormal clprobdistNormal;

/*! @name Functions for standard normal distribution
 *
 *  These functions are for the normal distribution with \f$\mu = 0\f$ and
 *  \f$\sigma = 1\f$, and they take neither a distribution object nor a list of
 *  distribution parameters.
 */
/*! @{ */

/*! @copydoc clprobdistDensity
 */
cl_double clprobdistStdNormalDensity(cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistCDF
 */
cl_double clprobdistStdNormalCDF(cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDF
 */
cl_double clprobdistStdNormalComplCDF(cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDF
 */
cl_double clprobdistStdNormalInverseCDF(cl_double u, clprobdistStatus* err);

/*! @} */


//***********************************************************************
// Constructor and Destructor
//***********************************************************************

/*! @name Functions to create and destroy distribution objects
 *
 * @{
 */

/*! @brief @copybrief clprobdistCreate()
 *
 *  Create a new normal distribution object.
 *  Since this function allocates memory for the new distribution object;
 *  clprobdistDestroy() must be called to release the allocated memory.
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     New distribution object.
 */
clprobdistNormal* clprobdistNormalCreate(cl_double mu, cl_double sigma, size_t* bufSize, clprobdistStatus* err);

/*! @copydoc clprobdistDestroy()
 */
clprobdistStatus clprobdistNormalDestroy(clprobdistNormal* dist);

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
cl_double clprobdistNormalDensityWithObject(const clprobdistNormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistCDFWithObject()
 */
cl_double clprobdistNormalCDFWithObject(const clprobdistNormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDFWithObject()
 */
cl_double clprobdistNormalComplCDFWithObject(const clprobdistNormal* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDFWithObject()
 */
cl_double clprobdistNormalInverseCDFWithObject(const clprobdistNormal* dist, cl_double u, clprobdistStatus* err);

/*! @copydoc clprobdistMeanWithObject()
 */
cl_double clprobdistNormalMeanWithObject(const clprobdistNormal* dist, clprobdistStatus* err);

/*! @copydoc clprobdistVarianceWithObject()
 */
cl_double clprobdistNormalVarianceWithObject(const clprobdistNormal* dist, clprobdistStatus* err);

/*! @copydoc clprobdistStdDeviationWithObject()
 */
cl_double clprobdistNormalStdDeviationWithObject(const clprobdistNormal* dist, clprobdistStatus* err);

/*! @brief Return the mean.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return Mean \f$\mu\f$.
 */
cl_double clprobdistNormalGetMu(const clprobdistNormal* dist, clprobdistStatus* err);

/*! @brief Return the standard deviation.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return Standard deviation \f$\sigma\f$.
 */
cl_double clprobdistNormalGetSigma(const clprobdistNormal* dist, clprobdistStatus* err);

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
 *  @see clprobdistNormalDensityWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistNormalDensity(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDF()
 *
 *  @see clprobdistNormalCDFWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistNormalCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDF()
 *
 *  @see clprobdistNormalComplCDFWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistNormalComplCDF(cl_double mu, cl_double sigma, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDF()
 *
 *  @see clprobdistNormalInverseCDFWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
cl_double clprobdistNormalInverseCDF(cl_double mu, cl_double sigma, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMean()
 *
 *  @see clprobdistNormalMeanWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return	Mean of the distribution.
 */
cl_double clprobdistNormalMean(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVariance()
 *
 *  @see clprobdistNormalVarianceWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistNormalVariance(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviation()
 *
 *  @see clprobdistNormalStdDeviationWithObject()
 *
 *  @param[in]  mu	Value of the mean \f$\mu\f$.
 *  @param[in]  sigma	Value of the standard deviation \f$\sigma\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistNormalStdDeviation(cl_double mu, cl_double sigma, clprobdistStatus* err);

/*! @} */

#ifdef __cplusplus
}
#endif

#endif
