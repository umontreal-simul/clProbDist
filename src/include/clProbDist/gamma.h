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
#ifndef CLPROBDIST_GAMMADIST_H
#define CLPROBDIST_GAMMADIST_H

/*! @file gamma.h
 *  @brief API of the gamma distribution
 *
 *  Implementation of \ref clProbDist_template.h for the gamma distribution,
 *  adapted from \cite iLEC08j .
 *
 *  For the sake of uniformity, all functions of the gamma distribution API
 *  that take explicit distribution parameters as their first arguments,
 *  instead of a distribution object, also include an integer `decprec`
 *  argument, documented in \ref distributions.
 *  Some of them actually use the argument value, others just ignore it.  This
 *  is indicated in the individual descriptions of these functions.
 */


#include "clProbDist/clProbDist.h"
#include "clProbDist/continuous.h"

/*! @brief Gamma distribution object [**device**]
 *
 *  A structure that represents a gamma distribution object.
 */
typedef struct _clprobdistGamma clprobdistGamma;

/***********************************
* Constructor & destructor
***********************************/

/*! @name Functions to create and destroy distribution objects
 *
 * @{
 */

/*! @brief @copybrief clprobdistCreate()
 *
 *  Create a new gamma distribution object.
 *  Since this function allocates memory for the new distribution object;
 *  clprobdistDestroy() must be called to release the allocated memory.
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec Value of \f$d\f$.
 *  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     New distribution object.
 */
clprobdistGamma * clprobdistGammaCreate(cl_double alpha, cl_double lambda, int decprec, size_t* bufSize, clprobdistStatus* err);

/*! @copydoc clprobdistDestroy()
 */
clprobdistStatus clprobdistGammaDestroy(clprobdistGamma* dist);

/*! @} */


// Dynamic functions

/*! @name Functions for use with explicit distribution parameters
 *
 *  @{
 */

/*! @copydoc clprobdistDensityWithObject()
 */
cl_double clprobdistGammaDensityWithObject(const clprobdistGamma* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistCDFWithObject()
 */
cl_double clprobdistGammaCDFWithObject(const clprobdistGamma* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDFWithObject()
 */
cl_double clprobdistGammaComplCDFWithObject(const clprobdistGamma* dist, cl_double x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDFWithObject()
 */
cl_double clprobdistGammaInverseCDFWithObject(const clprobdistGamma* dist, cl_double u, clprobdistStatus* err);

/*! @copydoc clprobdistMeanWithObject()
 */
cl_double clprobdistGammaMeanWithObject(const clprobdistGamma* dist, clprobdistStatus* err);

/*! @copydoc clprobdistVarianceWithObject()
 */
cl_double clprobdistGammaVarianceWithObject(const clprobdistGamma* dist, clprobdistStatus* err);

/*! @copydoc clprobdistStdDeviationWithObject()
 */
cl_double clprobdistGammaStdDeviationWithObject(const clprobdistGamma* dist, clprobdistStatus* err);

/*! @brief Return the shape parameter of the distribution  [**device**]
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return Value of the shape parameter \f$\alpha\f$.
 */
cl_double clprobdistGammaGetAlpha(const clprobdistGamma* dist, clprobdistStatus* err);

/*! @brief Return the scale parameter of the distribution  [**device**]
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return Value of the scale parameter \f$\lambda\f$.
 */
cl_double clprobdistGammaGetLambda(const clprobdistGamma* dist, clprobdistStatus* err);

/*! @} */




/*! @name Functions for use with explicit distribution parameters
 *
 *  @{
 */

/*! @brief @copybrief clprobdistDensity()
 *
 *  @see clprobdistGammaDensityWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec This parameter is unused but present for uniformity.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistGammaDensity(cl_double alpha, cl_double lambda, int decprec, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDF()
 *
 *  @see clprobdistGammaCDFWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec Value of \f$d\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistGammaCDF(cl_double alpha, cl_double lambda, int decprec, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDF()
 *
 *  @see clprobdistGammaComplCDFWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec Value of \f$d\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistGammaComplCDF(cl_double alpha, cl_double lambda, int decprec, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDF()
 *
 *  @see clprobdistGammaInverseCDFWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec Value of \f$d\f$.
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
cl_double clprobdistGammaInverseCDF(cl_double alpha, cl_double lambda, int decprec, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMean()
 *
 *  @see clprobdistGammaMeanWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec This parameter is unused but present for uniformity.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Mean of the distribution.
 */
cl_double clprobdistGammaMean(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVariance()
 *
 *  @see clprobdistGammaVarianceWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec This parameter is unused but present for uniformity.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistGammaVariance(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviation()
 *
 *  @see clprobdistGammaStdDeviationWithObject()
 *
 *  @param[in]  alpha   Value of the shape parameter \f$\alpha\f$.
 *  @param[in]  lambda  Value of the scale parameter \f$\lambda\f$.
 *  @param[in]  decprec This parameter is unused but present for uniformity.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistGammaStdDeviation(cl_double alpha, cl_double lambda, int decprec, clprobdistStatus* err);

/*! @} */
/*! @privatesection
 */

/*! @brief @copybrief clprobdistGammaCDF()
 *
 *  Same as clprobdistGammaCDF(), but for \f$\lambda = 1\f$.
 */
cl_double clprobdistGammaCDF_1(cl_double alpha, int d, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistGammaComplCDF()
 *
 *  Same as clprobdistGammaComplCDF(), but for \f$\lambda = 1\f$.
 */
cl_double clprobdistGammaComplCDF_1(cl_double alpha, int d, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistGammaInverseCDF()
 *
 *  Same as clprobdistGammaInverseCDF(), but for \f$\lambda = 1\f$.
 */
cl_double clprobdistGammaInverseCDF_1(cl_double alpha, int d, cl_double u, clprobdistStatus* err);

#endif
