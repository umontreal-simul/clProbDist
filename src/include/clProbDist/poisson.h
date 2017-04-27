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
#ifndef CLPROBDIST_POISSONDIST_H
#define CLPROBDIST_POISSONDIST_H

/*! @file poisson.h
 *  @brief API of the Poisson distribution
 *
 *  Implementation of \ref clProbDist_template.h for the Poisson distribution,
 *  adapted from \cite iLEC08j .
 */

#include "clProbDist/clProbDist.h"
#include <float.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*! @brief Poisson distribution object [**device**]
 *
 *  A structure that represents a Poisson distribution object.
 */
typedef struct _clprobdistPoisson clprobdistPoisson;

/***********************************
* Constructor & destructor
***********************************/


/*! @name Functions to create and destroy distribution objects
 *
 * @{
 */

/*! @brief @copybrief clprobdistCreate()
 *
 *  Create a new Poisson distribution object.
 *  Since this function allocates memory for the new distribution object;
 *  clprobdistDestroy() must be called to release the allocated memory.
 *
 *  @param[in] lambda  Value of the mean \f$\lambda\f$.
 *  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     New distribution object.
 */
clprobdistPoisson* clprobdistPoissonCreate(cl_double lambda, size_t* bufSize, clprobdistStatus* err);

/*! @copydoc clprobdistDestroy()
 */
clprobdistStatus clprobdistPoissonDestroy(clprobdistPoisson* dist);

/*! @} */

/***********************************
* Implementation
***********************************/


//Dynamic Functions

/*! @name Functions for use with a distribution object
 *
 * @{
 */

/*! @copydoc clprobdistProbWithObject()
 */
cl_double clprobdistPoissonProbWithObject(const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err);

/*! @copydoc clprobdistCDFWithObject()
 */
cl_double clprobdistPoissonCDFWithObject(const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err);

/*! @copydoc clprobdistComplCDFWithObject()
 *
 *  @warning The complementary distribution function is defined as
 *  \f$\bar F(j) = \mathbb P[X \geq j]\f$.
 */
cl_double clprobdistPoissonComplCDFWithObject(const clprobdistPoisson* dist, cl_int x, clprobdistStatus* err);

/*! @copydoc clprobdistInverseCDFWithObject()
 */
cl_int clprobdistPoissonInverseCDFWithObject(const clprobdistPoisson* dist, cl_double u, clprobdistStatus* err);

/*! @copydoc clprobdistMeanWithObject()
 */
cl_double clprobdistPoissonMeanWithObject(const clprobdistPoisson* dist, clprobdistStatus* err);

/*! @copydoc clprobdistVarianceWithObject()
 */
cl_double clprobdistPoissonVarianceWithObject(const clprobdistPoisson* dist, clprobdistStatus* err);

/*! @copydoc clprobdistStdDeviationWithObject()
 */
cl_double clprobdistPoissonStdDeviationWithObject(const clprobdistPoisson* dist, clprobdistStatus* err);

/*! @brief Return the value of the mean \f$\lambda\f$  [**device**]
 */
cl_double clprobdistPoissonGetLambda(const clprobdistPoisson* dist, clprobdistStatus* err);

/*! @} */


/*! @name Functions for use with explicit distribution parameters
 *
 *  @{
 */

/*! @brief @copybrief clprobdistProb()
 *
 *  @see clprobdistPoissonProbWithObject()
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistPoissonProb(cl_double lambda, cl_int x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDF()
 *
 *  @see clprobdistPoissonCDFWithObject()
 *
 *  @param[in]  lambda   Value of the mean \f$\lambda\f$.
 *  @param[in]  x        Value of \f$x\f$.
 *  @param[out] err	Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistPoissonCDF(cl_double lambda, cl_int x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDF()
 *
 *  @see clprobdistPoissonComplCDFWithObject()
 *
 *  @warning The complementary distribution function is defined as
 *  \f$\bar F(j) = \mathbb P[X \geq j]\f$.
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistPoissonComplCDF(cl_double lambda, cl_int x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDF()
 *
 *  @see clprobdistPoissonInverseCDFWithObject()
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
cl_int clprobdistPoissonInverseCDF(cl_double lambda, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMean()
 *
 *  @see clprobdistPoissonMeanWithObject()
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return	Mean of the distribution.
 */
cl_double clprobdistPoissonMean(cl_double lambda, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVariance()
 *
 *  @see clprobdistPoissonVarianceWithObject()
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistPoissonVariance(cl_double lambda, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviation()
 *
 *  @see clprobdistPoissonStdDeviationWithObject()
 *
 *  @param[in]  lambda  Value of the mean \f$\lambda\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistPoissonStdDeviation(cl_double lambda, clprobdistStatus* err);

/*! @} */


/*! @privatesection
 */

/*! @copydoc clprobdistInverseCDF()
 */
cl_int clprobdistDiscreteInverseCDF(const clprobdistPoisson* dist, cl_double u, clprobdistStatus* err);

#endif /* POISSONDIST_H */
