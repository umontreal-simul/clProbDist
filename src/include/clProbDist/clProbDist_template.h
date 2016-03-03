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
 *   David Munger <mungerd@iro.umontreal.ca>        (2015)
 *   Nabil Kemerchou <kemerchn@iro.umontreal.ca>    (2015)
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 */

#pragma once
#ifndef CLRNG_TEMPLATE_H
#define CLRNG_TEMPLATE_H

#error This file is a template for specific probability distributions; it should not be included as is.  Use a distribution-specific header such as poisson.h instead.

#include <clProbDist.h>


/*! @file clProbDist_template.h
 *  @brief Template of the API for specific probability distributions (not to be included as is!)
 *
 *  The function and type names in this API all start with \c clprobdist.
 *  In each specific implementation, prefixes to function names are expanded to
 *  a specific prefix and the generic distribution object type \c
 *  clprobdistObject is replaced with a specific distribution object type.
 *  Also, \ref vartype is replaced with `cl_double` for continuous
 *  distributions or with `cl_int` for discrete distributions.
 *  For example, the generic declaration
 *  \code
 *  vartype clprobdistInverseCDFWithObject(const clprobdistObject* dist, cl_double u, clprobdistStatus* err);
 *  \endcode
 *  expands to
 *  \code
 *  cl_int clprobdistPoissonInverseCDFWithObject(const clprobdistPoisson* dist, cl_double u, clprobdistStatus* err);
 *  \endcode
 *  for the Poisson distribution, and to
 *  \code
 *  cl_double clprobdistExponentialInverseCDFWithObject(const clprobdistExponential* dist, cl_double u, clprobdistStatus* err);
 *  \endcode
 *  for the normal distribution.
 *
 *  The API of each distribution, except for the standard normal that has no
 *  parameter, comes in two flavors: using a distribution object or using
 *  explicit distribution parameters.
 *  Since the number of parameters and their names vary across distributions,
 *  in the generic description of the API below, the list of distribution
 *  parameters is abbreviated by `DIST_PARAMS` in the function signatures, and
 *  should be replaced by the appropriate list of parameters, given in \ref
 *  distributions for each distribution.
 *
 *
 *  ### Host and Device APIs
 *
 *  The functions described here are all available on the host, in all implementations, 
 *  unless specified otherwise.  Only some of the functions and types are also
 *  available on the device in addition to the host;  they are tagged with
 *  [**device**].
 *  Other functions are only available on the device; they are tagged with
 *  [**device-only**].
 *  Some functions return an error code in \c err.
 *
 *  \todo Explain that on the device, the distribution object needs to be
 *  copied to local memory.  (Is this actually faster than the cache for global
 *  read-only memory?)
 */


/*! @brief Distribution object [**device**]
 *
 *  A structure that represents a distribution object.
 *  The contents of the structure depends on the type of distribution.
 *  It typically stores the parameters of the distribution, and possibly
 *  more information such as precomputed values.
 */
typedef struct { /* ... */ } clprobdistObject;

/*! @brief Random variable data type
 *
 *  Either \c cl_doubld or \c cl_int, depending on the distribution.
 *  See \ref distributions.
 *
 *  \todo Maybe remove this.  This is here to allow to document both continuous
 *  and discrete distributions together.
 */
typedef cl_double vartype;


/*! @name Functions to create and destroy distribution objects
 *
 * @{
 */

/*! @brief Create a distribution object.
 *
 *  Create a new distribution object.
 *  Since this function allocates memory for the new distribution object;
 *  clprobdistDestroy() must be called to release the allocated memory.
 *
 *  For each distribution type, the token \ref DIST_PARAMS expands to a
 *  list of the probability distribution parameters specific to the type of
 *  distribution.  See \ref distributions for more details.
 *
 *  @param[out] bufSize Size in bytes of the created distribution object, or \c NULL.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     New distribution object.
 */
clprobdistObject* clprobdistCreate(DIST_PARAMS, size_t* bufSize, clprobdistStatus* err);

/*! @brief Destroy a distribution object.
 *
 *  Release the resources associated to a distribution object.
 *
 *  @param[in,out] dist Distribution object.
 *  @return     Error status.
 */
clprobdistStatus clprobdistDestroy(clprobdistObject* dist);

/*! @} */


/*! @name Functions for use with a distribution object
 *
 *  These functions take a distribution object as their first argument.
 *  Some distribution objects contain precomputed values that can accelerate
 *  multiple evaluations of certain functions.
 *  When this is not needed, one can use the API below of functions with
 *  explicit distribution parameters, that do not require the prior creation of
 *  a distribution object.
 *
 * @{
 */

/*! @brief Probability density function [**device**]
 *
 *  Return \f$f(x)\f$, the value at \f$x=\f$`x` of the probability density function
 *  associated with the distribution object \c dist.
 *
 *  This function is defined only for continuous distributions (see \ref distributions).
 *
 *  @param[in]  dist    Distribution object.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistDensityWithObject(const clprobdistObject* dist, cl_double x, clprobdistStatus* err);

/*! @brief Probability mass function [**device**]
 *
 *  Return \f$p(x)\f$, the probability of \f$x\f$ associated with the
 *  distribution object \c dist.
 *
 *  This function is defined only for discrete distributions (see \ref distributions).
 *
 *  @param[in]  dist    Distribution object.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$p(x)\f$.
 */
cl_double clprobdistProbWithObject(const clprobdistObject* dist, cl_int x, clprobdistStatus* err);

/*! @brief Cumulative density function [**device**]
 *
 *  Return \f$F(x)\f$, the value at \f$x=\f$`x` of the distribution function
 *  associated with the distribution object \c dist.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistCDFWithObject(const clprobdistObject* dist, vartype x, clprobdistStatus* err);

/*! @brief Complementary CDF or reliability function [**device**]
 *
 *  Return \f$\bar F(x)\f$, the value of the complementary distribution
 *  function associated with the distribution object \c dist.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistComplCDFWithObject(const clprobdistObject* dist, vartype x, clprobdistStatus* err);

/*! @brief Inverse cumulative density function [**device**]
 *
 *  Return \f$F^{-1}(u)\f$, the value at \f$u=\f$`u` of the inverse
 *  distribution function associated with the distribution object \c dist.
 *  The type of the return value is `cl_int` for a discrete distribution of
 *  `cl_double` for a continuous distribution.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
vartype clprobdistInverseCDFWithObject(const clprobdistObject* dist, cl_double u, clprobdistStatus* err);

/*! @brief Mean of the distribution [**device**]
 *
 *  Return the mean of the distribution associated with the distribution object \c dist.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Mean of the distribution.
 */
cl_double clprobdistMeanWithObject(const clprobdistObject* dist, clprobdistStatus* err);

/*! @brief Variance of the distribution [**device**]
 *
 *  Return the variance of the distribution associated with the distribution object \c dist.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistVarianceWithObject(const clprobdistObject* dist, clprobdistStatus* err);

/*! @brief Standard deviation of the distribution [**device**]
 *
 *  Return the standard deviation of the distribution associated with the distribution object \c dist.
 *
 *  @param[in]  dist    Distribution object.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistStdDeviationWithObject(const clprobdistObject* dist, clprobdistStatus* err);

/*! @brief Copy a distribution object into already allocated memory [**device-only**]
 *
 *  Copy the distribution object @c srcDist located in global memory into the
 *  buffer @c destDist located in local or private memory.
 *  This function *does not* allocate memory for the distribution object, 
 *  it assumes that this has already been done.
 *
 *  When the distribution API is configured to use distribution objects stored
 *  in global memory (the default), there is no need to copy distribution
 *  objects across memory types, since the kernel already receives them in
 *  global memory (see @ref mem_types).
 *  The same applies to constant memory.
 *  In such cases, this function is not defined.

 *  @param[out] destDist	Destination buffer into which to copy (its
 *				            content will be overwritten).
 *                          The qualifier `_CLPROBDIST_<DIST>_OBJ_MEM`, where
 *                          `<DIST>` is the uppercase distribution name, is
 *                          replaced with the OpenCL keywords `__private` or
 *                          `__local`, respectively.
 *  @param[in]  srcDist		Distribution object to be copied.
 *
 *  @return     Error status
 *
 *  @warning This function is available only on the device.
 */
clprobdistStatus clprobdistCopyOverFromGlobal(_CLPROBDIST_<DIST>_OBJ_MEM clprobdistObject* destDist, __global const clprobdistObject* srcDist);

/*! @} */


/*! @name Functions for use with explicit distribution parameters
 *
 *  These functions do not require a previously created distribution objects.
 *  They take the distribution parameters as their first arguments
 *  (see the description of \ref DIST_PARAMS in \ref distributions).
 *  They cannot take advantage of precomputed values for the distributions.
 * 
 *  @{
 */

/*! @brief @copybrief clprobdistDensityWithObject()
 *
 *  Same as clprobdistDensityWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  This function is defined only for continuous distributions (see \ref distributions).
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistDensity(DIST_PARAMS, cl_double x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistProbWithObject()
 *
 *  Same as clprobdistProbWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  This function is defined only for discrete distributions (see \ref distributions).
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$f(x)\f$.
 */
cl_double clprobdistProb(DIST_PARAMS, cl_int x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistCDFWithObject()
 *
 *  Same as clprobdistCDFWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F(x)\f$.
 */
cl_double clprobdistCDF(DIST_PARAMS, vartype x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistComplCDFWithObject()
 *
 *  Same as clprobdistComplCDFWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[in]  x       Value of \f$x\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$\bar F(x)\f$.
 */
cl_double clprobdistComplCDF(DIST_PARAMS, vartype x, clprobdistStatus* err);

/*! @brief @copybrief clprobdistInverseCDFWithObject()
 *
 *  Same as clprobdistInverseCDFWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[in]  u       Value of \f$u \in [0,1]\f$.
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Value of \f$F^{-1}(u)\f$.
 */
vartype clprobdistInverseCDF(DIST_PARAMS, cl_double u, clprobdistStatus* err);

/*! @brief @copybrief clprobdistMeanWithObject()
 *
 *  Same as clprobdistMeanWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Mean of the distribution.
 */
cl_double clprobdistMean(DIST_PARAMS, clprobdistStatus* err);

/*! @brief @copybrief clprobdistVarianceWithObject()
 *
 *  Same as clprobdistVarianceWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Variance of the distribution.
 */
cl_double clprobdistVariance(DIST_PARAMS, clprobdistStatus* err);

/*! @brief @copybrief clprobdistStdDeviationWithObject()
 *
 *  Same as clprobdistStdDeviationWithObject(), but with explicit distribution
 *  parameters instead of distribution object.
 *
 *  See \ref distributions for the expansion of \ref DIST_PARAMS.
 *
 *  @param[out] err     Error status variable, or \c NULL.
 *  @return     Standard deviation of the distribution.
 */
cl_double clprobdistStdDeviation(DIST_PARAMS, clprobdistStatus* err);

/*! @} */

#endif /* CLPROBDIST_H *
* vim: syntax=c.doxygen spell spelllang=en fdm=syntax fdls=0 expandtab
*/
