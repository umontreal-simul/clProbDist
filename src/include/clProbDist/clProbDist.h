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

/*! @mainpage Introduction
 *
 *  We introduce clProbDist, an OpenCL library for probability distributions.
 *  It provides facilities for evaluating probability density or mass functions,
 *  distribution functions and their inverses, and reliability functions,
 *  currently for five distributions: normal, lognormal, exponential, gamma and Poisson.
 *  Applications include nonuniform variate generation, computing confidence
 *  intervals, computing \f$p\f$-values for statistical tests,
 *  computing densities to evaluate the likelihood function of a sample, etc.
 *
 *  This library, in its current state, defines a framework for working with
 *  probability distributions in OpenCL, and to show some working examples.
 *  Dozens of additional distributions need to be contributed to clProbDist to
 *  make it a mature library.
 *
 *  clProbDist allows users to generate nonuniform variates (see
 *  @ref example_poissongen for illustration) through inversion only.
 *  Other (possibly faster) methods, including acceptance-rejection based
 *  approaches, need to be implemented somewhere else.
 *  The inversion method relies on evaluating the inverse cumulative
 *  distribution function (CDF) of a given distribution.
 *  For certain distributions, it may be profitable to compute
 *  lookup arrays in advance and use them to evaluate the inverse CDF repeatedly
 *  with the same parameters, as when generating many random variates from the
 *  same distribution.
 *  clProbDist provides facilities to do so (see @ref example_poissongen_obj
 *  for an example).
 *
 *
 *  ### What to read next?
 *
 *  - @ref distributions
 *  - @ref basic_usage
 *  - @link clProbDist_template.h
 *    API supported by each probability distribution
 *    @endlink
 *    presented in a generic fashion
 *  - @link clProbDist.h
 *    API not related to specific probability distributions
 *    @endlink
 *  - @ref configuration
 *  - @ref mem_types (advanced usage)
 *
 *
 *
 *  @section basic_usage Basic usage
 *
 *  @subsection prefixes Distributions, prefixes and headers
 *
 *  The API is almost the same for every distribution.
 *  It is described in detail in clProbDist_template.h.
 *  Each distribution has its own header file (which is normally the lowercase
 *  name of the distribution, with a `.h` or `.clh` extension on the host or
 *  device, respectively), and each function contains the
 *  name of the distribution in its prefix (listed in @ref distributions).
 *  In the
 *  @link clProbDist_template.h generic probability distribution API reference@endlink,
 *  the distribution-specific part of the prefix is not shown.
 *
 *  For example, functions from the @ref normal_dist API have names that begin
 *  with `clprobdistNormal` and are declared in the normal.h header file.
 *  In particular, the probability density function clprobdistDensity(), the
 *  cumulative distribution function clprobdistCDF() and its inverse
 *  clprobdistInverseCDF() from the generic API, expand to
 *  clprobdistNormalDensity(), clprobdistNormalCDF() and
 *  clprobdistNormalInverseCDF(), respectively, for the normal distribution.
 *  More concretely, to store in `p` the probability density of `x` in a normal
 *  distribution of mean \f$\mu\f$ = `mu` and standard deviation
 *  \f$\sigma\f$ = `sigma`, one would invoke
 *  @code
 *  cl_double p = clprobdistNormalDensity(mu, sigma, x, NULL);
 *  @endcode
 *  Note that `x`, `mu` and `sigma` in the above are all floating-point values
 *  of type `cl_double`.
 *  The last argument can be used for reporting errors, but errors are ignored
 *  here by passing `NULL`.
 *
 *
 *  @subsection example_ci Example 1: computing a confidence interval
 *
 *  The following code defines a function (on the host) that computes a normal
 *  confidence interval of level `level` on the mean of a set of
 *  `n` observations of type `cl_double` stored in the array `observations`.
 *  This assumes, of course, that the mean is normally distributed.
 *  The lower and upper bounds are returned in the `lower` and `upper`
 *  variables, respectively:
 *  @code
 *  #include <clProbDist/normal.h>
 *  #include <math.h> // for sqrt()
 *
 *  void normalConfidenceInterval(cl_int n, cl_double* observations, cl_double level,
 *                                cl_double* lower, cl_double* upper)
 *  {
 *      // compute sample average and variance
 *      cl_double sum = 0.0;
 *      cl_double sum_squares = 0.0;
 *      for (size_t i = 0; i < n; i++) {
 *          sum         += observations[i];
 *          sum_squares += observations[i] * observations[i];
 *      }
 *      cl_double average  = sum / n;
 *      cl_double variance = (sum_squares - average * sum) / (n - 1);
 *
 *      // compute normal distribution parameters
 *      cl_double mu       = average;
 *      cl_double sigma    = sqrt(variance / n);
 *      // compute confidence interval
 *      *lower = clprobdistNormalInverseCDF(mu, sigma, 0.5 * (1.0 - level), NULL);
 *      *upper = clprobdistNormalInverseCDF(mu, sigma, 0.5 * (1.0 + level), NULL);
 *  }
 *  @endcode
 *  @footnote{Because there will eventually be many distributions with many
 *  header files, I deemed it best to group them into a directory (named
 *  clProbDist).  On Linux systems, this helps to avoid cluttering the standard
 *  include directory.
 *  Numerous multi-header C libraries (OpenCL, Cairo, GNU Scientific Library, GTK+
 *  3.0, libxml 2, SDL, OpenGL) put their header files in a subdirectory named
 *  after the library.}
 *  The first part of `normalConfidenceInterval()` simply computes the sample average and
 *  variance of the observations.
 *  The second part illustrates how to evaluate the inverse normal distribution
 *  function by invoking clprobdistNormalInverseCDF() with the normal
 *  distribution parameters \f$\mu\f$=`mu` and \f$\sigma\f$=`sigma` (the mean
 *  and variance) as its first arguments, for two different quantiles.
 *  The last argument is used for reporting errors, which are ignored here by
 *  setting it to `NULL`.
 *
 *  A more efficient implementation would compute the half width of the
 *  interval using a single call to clprobdistStdNormalInverseCDF() from the
 *  @ref stdnormal_dist API.
 *  The standard normal distribution always has zero mean and a unit variance,
 *  so functions from its API do not take any distribution parameters.
 *  With this, we could replace the last four instructions of the
 *  normalConfidenceInterval defined above with:
 *  @code
 *      cl_double half_width = sqrt(variance / n) *
 *                clprobdistStdNormalInverseCDF(0.5 * (1.0 + level), NULL);
 *      *lower = average - half_width;
 *      *upper = average + half_width;
 *  @endcode
 *
 *  The complete code for this example can be found in @ref DocsTutorial/example1.c.
 *
 *
 *  @subsection example_poissongen Example 2: non-uniform variate generation
 *
 *  In this example, we want to simulate a stochastic system (not completely
 *  specified as our focus here is the usage of clProbDist) that uses a
 *  number, say 100, of Poisson variates, all with the same mean `lambda`, as
 *  its input.
 *  To generate a Poisson variate with the inversion method, we first generate
 *  a pseudorandom number uniformly distributed in (0,1) using the MRG31k3p
 *  generator from the
 *  [clRNG library](https://github.com/clMathLibraries/clRNG),
 *  and we map the output through clprobdistPoissonInverseCDF().
 *  (The user is referred to the
 *  [clRNG documentation](http://clmathlibraries.github.io/clRNG/htmldocs/index.html)
 *  and
 *  [clRNG tutorial document](http://clmathlibraries.github.io/clRNG/docs/clrng-api.pdf)
 *  for instructions on random streams usage.)
 *  For simplicity, we assume that the Poisson variates are used sequentially,
 *  one at a time, so we do not need to store them all at once in advance.
 *  Again, we run everything on the host.
 *  The code to average 1024 realizations of the output of this stochastic
 *  system would look like this:
 *  @code
 *  #include <clProbDist/poisson.h>
 *  #include <clRNG/mrg31k3p.h>
 *
 *  cl_double simulateOneRun(cl_double lambda, clrngMrg31k3pStream* stream)
 *  {
 *      cl_double output = 0.0;
 *      for (int i = 0; i < 100; i++) {
 *          cl_double u = clrngMrg31k3pRandomU01(stream);
 *          cl_int poisson_variate = clprobdistPoissonInverseCDF(lambda, u, NULL);
 *          // modify output using poisson_variate
 *          output += poisson_variate;
 *      }
 *      return output;
 *  }
 *
 *  int main()
 *  {
 *      // prepare
 *      int replications = 1024;
 *      cl_double lambda = 50.0;
 *      clrngMrg31k3pStream* stream = clrngMrg31k3pCreateStreams(NULL, 1, NULL, NULL);
 *
 *      // simulate
 *      cl_double sum = 0.0;
 *      for (int i = 0; i < replications; i++)
 *          sum += simulateOneRun(lambda, stream);
 *      printf("The average output is %.3f\n", sum / replications);
 *
 *      // clean up
 *      clrngMrg31k3pDestroyStreams(stream);
 *
 *      return 0;
 *  }
 *  @endcode
 *  Here, the output of `simulateOneRun()` is simply the sum of the 100 Poisson
 *  variates.
 *  In practice, it would be more elaborate, but the purpose of this example is
 *  to show how to generate nonuniform variates.
 *  The call to `clprobdistPoissonInverseCDF(lambda, u, NULL)` evaluates the
 *  inverse Poisson distribution function with mean `lambda` at quantile `u`,
 *  ignoring errors.
 *
 *  The complete code for this example can be found in @ref DocsTutorial/example2.c.
 *
 *  The above code generates 1024 × 100 Poisson variates with the same mean,
 *  and can be very inefficient when `lambda` is large.
 *  A faster approach makes use of a lookup array to avoid repeating costly computations.
 *  Such facility is provided by a selection of functions of the Poisson
 *  distribution API, illustrated in @ref example_poissongen_obj below.
 *  For this purpose, we introduce @ref dist_objects.
 *
 *
 *  @subsection dist_objects Distribution objects
 *
 *  The API of each distribution comes in two variants.
 *  One uses distribution parameters as the first arguments of the functions;
 *  the other replaces the distribution parameters with a
 *  distribution object, and appends `WithObject` to the function names.
 *  @footnote{Alternatively, functions that take a distribution object could be
 *  named without the `WithObject` suffix whereas functions that take distribution
 *  parameters could be appended `WithParams` to their names.  But, this would
 *  be inappropriate for the API of the standard normal distribution that takes
 *  no explicit parameters (the mean and standard deviation are implicitly set
 *  to 0 and 1, respectively).}
 *
 *  The type of a distribution object is specific to the distribution it
 *  represents, e.g., the Poisson distribution has objects of type
 *  `clprobdistPoisson`.
 *  A distribution object stores the parameters of the distribution it
 *  represents, and, for some distributions, it may also store additional data.
 *
 *  In particular, for discrete distributions such as the Poisson one,
 *  the distribution object also stores lookup arrays to speed up the
 *  evaluation of the probability, of the distribution function and its
 *  inverse, and of the reliability function.
 *  Building the lookup arrays cost some time, so it is advisable to use a
 *  distribution object especially when these functions need to be evaluated
 *  many times.
 *  Before before using a distribution object with these functions, it must be
 *  created using clprobdistCreate(); it can be destroyed with clprobdistDestroy().
 *
 *  For other distributions, it should make very little difference in terms of
 *  performance which variant of the distribution API is used (with
 *  distribution objects or distribution parameters).
 *
 *  In the
 *  @link clProbDist_template.h generic probability distribution API reference@endlink,
 *  the generic type name for a distribution object is clprobdistObject.
 *
 *
 *  @subsection example_poissongen_obj Example 3: non-uniform variate generation revisited
 *
 *  Here, we modify @ref example_poissongen to use a distribution object in
 *  order to accelerate the multiple evaluations of the inverse distribution
 *  function.
 *  First, `simulateOneRun()` needs to receive a distribution object instead of
 *  the mean of the Poisson distribution, so its signature becomes
 *  @code
 *  cl_double simulateOneRun(const clprobdistPoisson* dist, clrngMrg31k3pStream* stream)
 *  @endcode
 *  Next, we replace the call to clprobdistPoissonInverseCDF() to its variant
 *  that takes a distribution object:
 *  @code
 *  cl_int poisson_variate = clprobdistPoissonInverseCDFWithObject(dist, u, NULL);
 *  @endcode
 *  Finally, we create the distribution object with:
 *  @code
 *  clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, NULL, NULL);
 *  @endcode
 *  and pass it to `simulateOneRun()` instead of the value of the mean:
 *  @code
 *  sum += simulateOneRun(dist, stream);
 *  @endcode
 *  The first argument to clprobdistPoissonCreate() is the value of distribution
 *  parameter \f$\lambda\f$.
 *  For distributions with multiple parameters, their values would be listed as
 *  the first arguments to clprobdistCreate() (expanded to the proper name for
 *  the distribution).
 *  For instance, to create a distribution object for the normal distribution
 *  with mean `mu` and standard deviation `sigma`, we would call:
 *  @code
 *  clprobdistNormal* dist = clprobdistNormalCreate(mu, sigma, NULL, NULL);
 *  @endcode
 *  The last two arguments to clprobdistPoissonCreate(), or to
 *  clprobdistCreate() in general, can be used to return the size of
 *  the created distribution object and the error status.
 *  Here, we ignore them by setting them to `NULL`.
 *
 *  The complete code for this example can be found in @ref DocsTutorial/example3.c.
 *
 *
 *  @subsection device_usage Usage from device code
 *
 *  Using clProbDist from device code is very similar to using it from host
 *  code.
 *  The device API variant that uses distribution parameters is the same as the
 *  host API.
 *  For example, to store in `p` the probability density of `x` in a normal
 *  distribution of mean \f$\mu\f$ = `mu` and standard deviation
 *  \f$\sigma\f$ = `sigma`, one would include the header file normal.clh and
 *  invoke:
 *  @code
 *  double p = clprobdistNormalDensity(mu, sigma, x);
 *  @endcode
 *  Here, `x`, `mu` and `sigma` are all floating-point values of type `double`.
 *
 *  The device API has no functions to create or destroy distribution objects.
 *  Besides that, it is almost the same as the host API.
 *  Here is an example kernel that receives a normal distribution object and
 *  uses it to invoke clprobdistNormalDensity():
 *  @code
 *  __kernel void my_kernel(__global const clprobdistNormal* dist)
 *  {
 *      double p = clprobdistNormalDensity(dist, 8.0, NULL);
 *      ...
 *  }
 *  @endcode
 *  Again, the last argument to clprobdistNormalDensity(), that can be used for
 *  error reporting but is ignored here by setting it to `NULL`.
 *  It is possible to store the distribution objects into other memory types
 *  than global memory, but this requires additional instructions in the host
 *  and kernel codes.
 *  This is discussed in @ref mem_types.
 *
 *
 *  @subsection example_poissongen_device Example 4: Poisson variate generation on the device
 *
 *  Here, we adapt @ref example_poissongen_obj, which runs on the host, to run
 *  on the device.
 *  First, we translate the computation code into device code, by replacing the
 *  extensions of the header names with `.clh` and by replacing `cl_double`
 *  with `double`, and by writing a kernel to receive data from the host, call
 *  `simulateOneRun()` and send data back to the host:
 *  @code
 *  #include <clProbDist/poisson.clh>
 *  #include <clRNG/mrg31k3p.clh>
 *
 *  double simulateOneRun(__global const clprobdistPoisson* dist, clrngMrg31k3pStream* stream)
 *  {
 *      double output = 0.0;
 *      for (int i = 0; i < 100; i++) {
 *          double u = clrngMrg31k3pRandomU01(stream);
 *          int poisson_variate = clprobdistPoissonInverseCDFWithObject(dist, u, NULL);
 *          // modify output using poisson_variate
 *          ...
 *      }
 *      return output;
 *  }
 *
 *  __kernel void my_kernel(
 *  	__global const clprobdistPoisson* dist,
 *  	__global const clrngMrg31k3pHostStream* host_streams,
 *  	__global double* output)
 *  {
 *      size_t gid = get_global_id(0);
 *      clrngMrg31k3pStream stream;
 *      clrngMrg31k3pCopyOverStreamsFromGlobal(1, &stream, &host_streams[gid]);
 *      output[gid] = simulateOneRun(dist, &stream);
 *  }
 *  @endcode
 *  The `__global` OpenCL C keyword in front of `const clprobdistPoisson*`
 *  indicates that the Poisson distribution object is stored in the device's
 *  global memory.
 *
 *  The host code is responsible for creating the distribution object and
 *  sending it into the device's memory.
 *  The latter requires the size of the object to be known, so we allow
 *  clprobdistPoissonCreate() to return it in its second argument (`dist_buf_size`):
 *  @code
 *  size_t dist_buf_size;
 *  clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, &dist_buf_size, NULL);
 *  @endcode
 *  After this instruction, `dist_buf_size` contains the size in bytes of the
 *  distribution object pointed to by `dist`.
 *  The size of the array is determined at runtime; it may vary depending on
 *  the value of the distribution parameter \f$\lambda\f$, set to `50.0` here.
 *  An OpenCL buffer must be created to copy the distribution object into the
 *  device's memory:
 *  @code
 *  cl_mem dist_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
 *                                   dist_buf_size, dist, NULL);
 *  @endcode
 *  In the above code snippet, the `context` variable is an OpenCL context.
 *  Finally, the buffer must be set as the first argument of the kernel (stored
 *  in the `kernel` variable):
 *  @code
 *  clSetKernelArg(kernel, 0, sizeof(dist_buf), &dist_buf);
 *  @endcode
 *
 *  The complete code for this example can be found in @ref DocsTutorial/example4.c
 *  and @ref DocsTutorial/example4_kernel.cl.
 *
 *
 *  @section configuration Configuration
 *
 *  @subsection environment Environment variables
 *
 *  For all features of the library to work properly, the `CLPROBDIST_ROOT`
 *  environment variable must be set to point to the installation path of the
 *  clProbDist package, that is, the directory under which lies the
 *  `include/clProbDist` subdirectory.
 *  Means of setting an environment variable depend on the operating system
 *  used.
 *
 *
 *  @section mem_types Device memory types
 *
 *  The API of each distribution assumes that the distribution objects are
 *  stored in a specific type of memory.
 *  It defaults to global memory, but can be customized by the user by changing
 *  the value of the preprocessor symbol `CLPROBDIST_<DIST>_OBJ_MEM`, where
 *  `<DIST>` is the uppercase name of the distribution, before including the
 *  device header file of the distribution, to one of the following values:
 *
 *  - CLPROBDIST_MEM_TYPE_PRIVATE to use private memory (except for the Poisson distribution);
 *  - CLPROBDIST_MEM_TYPE_LOCAL to use local memory;
 *  - CLPROBDIST_MEM_TYPE_CONSTANT to use constant memory;
 *  - CLPROBDIST_MEM_TYPE_GLOBAL to use global memory (the **default**).
 *
 *  For example, to store exponential distribution objects in **constant memory**,
 *  the device code should simply begin with:
 *  @code
 *  #define CLPROBDIST_EXPONENTIAL_OBJ_MEM CLPROBDIST_MEM_TYPE_CONSTANT
 *  #include <clProbDist/exponential.clh>
 *  @endcode
 *
 *  @anchor clprobdistCopyOverFromGlobal
 *  For other memory types than global and constant memory, the device API
 *  provides a `clprobdistCopyOverFromGlobal()` function, where the
 *  `clprobdist` prefix has to be expanded with the distribution name as
 *  explained in @ref prefixes.
 *  We give examples below of how to use this function.
 *
 *  To store exponential distribution objects in **private memory** instead,
 *  in addition to replacing `CLPROBDIST_MEM_TYPE_CONSTANT` with
 *  `CLPROBDIST_MEM_TYPE_PRIVATE`, each work item needs to make a private copy
 *  of the global distribution object and use the address of that copy with API
 *  calls, e.g.,
 *  @code
 *  #define CLPROBDIST_EXPONENTIAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE
 *  #include <clProbDist/exponential.clh>
 *  __kernel void example(__global const clprobdistExponential* gdist)
 *  {
 *      // allocate private storage for distribution object
 *      clprobdistExponential pdist;
 *      // copy distribution object from global memory to private memory
 *      clprobdistExponentialCopyOverFromGlobal(&pdist, gdist);
 *      // following API calls should use &pdist instead of gdist
 *      ...
 *  }
 *  @endcode
 *
 *  To store distribution objects in **local memory**, the user must
 *  first reserve local memory from the host code.
 *  For example, with a Poisson distribution object, the host code would
 *  include something like the following:
 *  @code
 *  // create a Poisson distribution object
 *  clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, NULL, NULL);
 *  // create OpenCL buffer for distribution object
 *  size_t dist_size;
 *  dist_buf = clCreateBuffer(context,
 *    CL_MEM_HOST_WRITE_ONLY | CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
 *    dist_size, dist, &err);
 *  // send the distribution object in global memory
 *  clSetKernelArg(kernel, 0, sizeof(dist_buf), &dist_buf);
 *  // allocate local memory to allow the device to copy the distribution object
 *  clSetKernelArg(kernel, 1, dist_size, NULL);
 *  @endcode
 *  In the kernel, the distribution object must be copied from global memory to
 *  local memory, which is shared by all work items from the same work group.
 *  Perhaps the simplest way to do so is to have the first work item of each
 *  work group do the copy, e.g.,
 *  @code
 *  #include <clProbDist/clProbDist.clh>
 *  #define CLPROBDIST_POISSON_OBJ_MEM CLPROBDIST_MEM_TYPE_LOCAL
 *  #include <clProbDist/poisson.clh>
 *  __kernel void example(__global const clprobdistPoisson* gdist,
 *                        __local        clprobdistPoisson* ldist)
 *  {
 *      int gid = get_global_id(0);
 *      if (get_global_id(0) % get_local_size(0) == 0)
 *          // copy distribution object from global memory to local memory
 *          clprobdistPoissonCopyOverFromGlobal(ldist, gdist);
 *      // wait for all work items of the work group to be ready
 *      barrier(CLK_LOCAL_MEM_FENCE);
 *      // following API calls should use ldist instead of gdist
 *      ...
 *  }
 *  @endcode
 *  This is substantially more complicated than using global memory.  It may be
 *  faster or slower depending on the distribution and its parameters, the
 *  OpenCL hardware and the vendor implementation of OpenCL.
 *
 *
 *  @section distributions Implemented distributions
 *
 *  The following table lists the distributions that are currently implemented
 *  in clProbDist with the type of their support (continuous or discrete), the prefix
 *  for functions of the distribution API, and the name of the corresponding
 *  header file on the host and on the device.
 *
 *  | Distribution    | Support Type  | Prefix                  | Host Header File | Device Header File |
 *  | --------------  | ------------- | ----------------------- | ---------------- | ------------------ |
 *  | Normal          | continuous    | `clprobdistNormal`      | normal.h         | normal.clh         |
 *  | Standard Normal | continuous    | `clprobdistStdNormal`   | normal.h         | normal.clh         |
 *  | Lognormal       | continuous    | `clprobdistLognormal`   | lognormal.h      | lognormal.clh      |
 *  | Exponential     | continuous    | `clprobdistExponential` | exponential.h    | exponential.clh    |
 *  | Gamma           | continuous    | `clprobdistGamma`       | gamma.h          | gamma.clh          |
 *  | Poisson         | discrete      | `clprobdistPoisson`     | poisson.h        | poisson.clh        |
 *
 *  The API of continuous distributions provide a clprobdistDensity()
 *  function to compute the probability density; for discrete distributions, it
 *  is replaced by a clprobdistProb() function that computes the probability
 *  mass.
 *
 *  The common API for all distributions, in a generic form, is given in @ref
 *  clProbDist_template.h.
 *  Some parts of the API vary across distributions, for instance, the type of
 *  return value of the inverse distribution function: `cl_double` for continuous
 *  distributions and `cl_int` for discrete distributions.
 *  This is in fact the type of the random variable described by the
 *  distribution, and its programmatic type is denoted by `vartype` in the
 *  generic API.
 *  Also, the number of distribution parameters vary, as well as their
 *  individual types and names.
 *  The token `DIST_PARAMS` substitutes for the list of distribution parameters
 *  in the generic API.
 *  The table below explicits the expansion of `vartype` and `DIST_PARAMS` for
 *  implemented distributions.
 *  The meaning of the distribution parameters is given below for each
 *  individual distribution.
 *
 *  | Distribution    | `vartype`   | `DIST_PARAMS` @anchor DIST_PARAMS                                     |
 *  | --------------  | ----------- | --------------------------------------------------------------------- |
 *  | Normal          | `cl_double` | `cl_double` \f$\mu\f$,    `cl_double` \f$\sigma\f$                    |
 *  | Standard Normal | `cl_double` | *none*                                                                |
 *  | Lognormal       | `cl_double` | `cl_double` \f$\mu\f$,    `cl_double` \f$\sigma\f$                    |
 *  | Exponential     | `cl_double` | `cl_double` \f$\lambda\f$ `cl_double`                                 |
 *  | Gamma           | `cl_double` | `cl_double` \f$\alpha\f$, `cl_double` \f$\lambda\f$, `cl_int` \f$d\f$ |
 *  | Poisson         | `cl_int`    | `cl_double` \f$\lambda\f$                                             |
 *
 *  In addition to the functions listed with generic API names in
 *  clProbDist_template.h, each parameterized distribution exposes functions to
 *  retrieve the values of the parameters a distribution object was created
 *  with.
 *
 *
 *  @subsection normal_dist Normal distribution
 *
 *  With mean \f$\mu\f$ and variance \f$\sigma^2\f$, where \f$\sigma > 0\f$,
 *  its density function is
 *  \f[
 *    f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left[ \frac{-(x-\mu)^2}{2\sigma^2} \right]
 *    \qquad
 *    (-\infty < x < \infty).
 *  \f]
 *
 *  For the normal distribution, @ref DIST_PARAMS expands to:
 *  @code
 *  cl_double mu, cl_double sigma
 *  @endcode
 *  where `mu` \f$=\mu\f$ and `sigma` \f$=\sigma\f$.
 *  The API of the normal distribution adds the functions
 *  clprobdistNormalGetMu() and clprobdistNormalGetSigma() to retrieve the
 *  value of these parameters.
 *
 *  The CDF is evaluated using Uses the Chebyshev approximation proposed in
 *  \cite tSCH78a, which gives 16 decimals of precision.
 *  The inverse CDF is computed using different rational Chebyshev
 *  approximations \cite tBLA76a, which give 16 decimal digits of precision for
 *  \f$2.2 \times 10^{-308} < u < 1\f$.
 *
 *  The complete API of the normal distribution is described in normal.h.
 *
 *
 *  @subsection stdnormal_dist Standard normal distribution
 *
 *  This special case of the normal distribution with \f$\mu = 0\f$ and
 *  \f$\sigma = 1\f$ has distribution function
 *  \f[
 *    F(x) = \Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\infty}^x e^{-t^2/2} \,dt
 *    \qquad
 *    (-\infty < x < \infty).
 *  \f]
 *
 *  For the standard normal distribution, there are no free parameters, so
 *  @ref DIST_PARAMS expands to nothing.
 *
 *  The complete API of the standard normal distribution is also described in
 *  normal.h.
 *
 *
 *  @subsection lognormal_dist Lognormal distribution
 *
 *  With scale parameter \f$\mu\f$ and shape parameter \f$\sigma^2\f$, where \f$\sigma > 0\f$,
 *  its density function is
 *  \f[
 *    f(x) = \frac{1}{\sqrt{2 \pi} \sigma x} \exp\left[ \frac{-(\ln x - \mu)^2}{2\sigma^2} \right]
 *    \qquad
 *    (0 < x < \infty).
 *  \f]
 *  Its distribution function is
 *  \f[
 *    F(x) = \Phi((\ln x - \mu) / \sigma)
 *    \qquad
 *    (0 < x < \infty),
 *  \f]
 *  where \f$\Phi\f$ is the distribution function for the @ref stdnormal_dist.
 *  This inverse distribution function is given by:
 *  \f[
 *    F^{-1}(u) = \exp\left[ \mu + \sigma\Phi^{-1}(u) \right]
 *    \qquad
 *    (0 \leq u < 1).
 *  \f]
 *  If \f$\ln Y\f$ has a normal distribution, then \f$Y\f$ has a lognormal
 *  distribution with the same parameters.
 *
 *  The CDF and inverse CDF are evaluated using the implementations for the
 *  normal distribution.
 *
 *  For the lognormal distribution, @ref DIST_PARAMS expands to:
 *  @code
 *  cl_double mu, cl_double sigma
 *  @endcode
 *  where `mu` \f$=\mu\f$ and `sigma` \f$=\sigma\f$.
 *  The API of the normal distribution adds the functions
 *  clprobdistLognormalGetMu() and clprobdistLognormalGetSigma() to retrieve the
 *  value of these parameters.
 *
 *  The complete API of the normal distribution is described in lognormal.h.
 *
 *
 *  @subsection exponential_dist Exponential distribution
 *
 *  With mean \f$1/\lambda\f$, where \f$\lambda > 0\f$, its density,
 *  distribution and inverse distribution functions are, respectively,
 *  \f{align}{
 *    f(x)      &= \lambda e^{-\lambda x}    & (x \geq 0) \\
 *    F(x)      &= 1 - e^{-\lambda x}        & (x \geq 0) \\
 *    F^{-1}(x) &= -\lambda^{-1} \ln(1 - u)  & (0 \leq u < 1).
 *  \f}
 *  The inverse CDF is computed exactly with the above formula.
 *
 *  For the exponential distribution, @ref DIST_PARAMS expands to:
 *  @code
 *  cl_double lambda
 *  @endcode
 *  where `lambda` \f$=\lambda\f$.
 *  The API of the exponential distribution adds the function
 *  clprobdistExponentialGetLambda() retrieve the value of this parameter.
 *
 *  The complete API of the exponential distribution is described in
 *  exponential.h.
 *
 *
 *  @subsection gamma_dist Gamma distribution
 *
 *  With shape parameter \f$\alpha > 0\f$ and scale parameter \f$\lambda >
 *  0\f$, its density function is
 *  \f[
 *    f(x) = \frac{\lambda^\alpha x^{\alpha - 1} e^{\lambda x}}{\Gamma(\alpha)}
 *    \qquad
 *    (x > 0),
 *  \f]
 *  where
 *  \f[
 *    \Gamma(\alpha) = \int_0^\infty x^{\alpha - 1} e^{-x} \,dx
 *  \f]
 *  is the usual gamma function.
 *  In particular, \f$\Gamma(n) = (n - 1)!\f$ when \f$n\f$ is a positive
 *  integer.
 *
 *  Several functions of the API for the gamma distribution are evaluated using
 *  an approximation with a precision of roughly \f$d\f$ decimal digits, where
 *  \f$d\f$ can be specified by the user.
 *  For consistency, \f$d\f$ must always be specified, but is ignored by
 *  clprobdistGammaDensity(), clprobdistGammaMean(),
 *  clprobdistGammaStdDeviation() and clprobdistGammaVariance().
 *
 *  The CDF is approximated with an improved version of the algorithm in
 *  \cite tBHA70a .
 *  The functions clprobdistGammaCDF() and clprobdistGammaCDFWithObject() try
 *  to return \f$d\f$ decimals digits of precision.  For \f$\alpha\f$ not
 *  too large (e.g., \f$\alpha \leq 1000\f$), \f$d\f$ gives a good idea of the
 *  precision attained.
 *  The inverse CDF \f$F^{-1}(u)\f$ is approximated by solving \f$F(x) - u =
 *  0\f$ for \f$x\f$ with the bisection method when \f$u \leq 10^{-8}\f$ or
 *  \f$\alpha \leq 1.5\f$, or with the Brent-Dekker method (see \cite iBRE71a
 *  and \cite iBRE73a) otherwise.
 *
 *  For the gamma distribution, @ref DIST_PARAMS expands to:
 *  @code
 *  cl_double alpha, cl_double lambda, cl_int decprec
 *  @endcode
 *  where `alpha` \f$=\alpha\f$, `lambda` \f$=\lambda\f$ and `decprec` \f$=d\f$.
 *  The API of the gamma distribution adds the functions
 *  clprobdistGammaGetAlpha() and clprobdistGammaGetLambda() to retrieve the
 *  value of the first two parameters.
 *
 *  The complete API of the gamma distribution is described in gamma.h.
 *
 *
 *  @subsection poisson_dist Poisson distribution
 *
 *  With mean \f$\lambda \geq 0\f$, its mass and distribution functions are,
 *  respectively,
 *  \f{align}{
 *    p(x) &= \frac{e^{-\lambda} \lambda^x}{x!}               & (x = 0, 1, \ldots) \\
 *    F(x) &= e^{-\lambda} \sum_{j=0}^x \frac{\lambda^j}{j!}  & (x = 0, 1, \ldots)
 *  \f}
 *
 *  When \f$\lambda > 200\f$, to compute \f$F(x)\f$ and \f$\bar F(x)\f$,
 *  clprobdistPoissonCDF() and clprobdistPoissonComplCDF() exploit the
 *  relationship \f$F(x) = 1 - G_{x+1}(\lambda)\f$, where \f$G_{x+1}\f$ is the
 *  gamma distribution function with parameters
 *  \f$(\alpha, \lambda) = (x + 1, 1)\f$.
 *  Otherwise, the probabilities are summed term by term.
 *  To evaluate the inverse CDF \f$F^{-1}(u)\f$, clprobdistPoissonInverseCDF()
 *  performs a sequential search for the smallest value of \f$x\f$ such that
 *  \f$F(x) \geq u\f$, by incrementing \f$x\f$ starting from \f$x = 0\f$.
 *  When \f$\lambda\f$ is very large, the probabilities
 *  \f$p(x)\f$ are empirically negligible when \f$x\f$ is far from
 *  \f$\lambda\f$.  So, when \f$\lambda \geq 700\f$, the implementation of
 *  clprobdistPoissonInverseCDF() starts with a binary search for the lower
 *  boundary \f$x_{\mathrm{min}} \in \{0,\dots,\lambda\}\f$ of an interval that
 *  contains the non-negligible probability values.  Then, the sequential
 *  search begins at \f$x = x_{\mathrm{min}}\f$.
 *
 *  If one has to compute \f$p(x)\f$ or \f$F(x)\f$ for several values of
 *  \f$x\f$ with the same \f$\lambda\f$, where \f$\lambda\f$ is not too large,
 *  then it is more efficient to instantiate an object and use the API variant
 *  that uses distribution objects of type clprobdistPoisson, e.g.,
 *  clprobdistPoissonProbWithObject() and clprobdistPoissonCDFWithObject(),
 *  for which return values are computed in advance and stored in lookup
 *  arrays.
 *  Probabilities smaller than \f$10^{-16}\f$ are not stored in the Poisson
 *  distribution object, but are computed directly each time they are needed
 *  (which should be very seldom).
 *  The inverse CDF is evaluated by clprobdistPoissonInverseCDFWithObject()
 *  with a binary search.
 *
 *  For the Poisson distribution, @ref DIST_PARAMS expands to:
 *  @code
 *  cl_double lambda
 *  @endcode
 *  where `lambda` \f$=\lambda\f$.
 *  The API of the exponential distribution adds the function
 *  clprobdistPoissonGetLambda() retrieve the value of this parameter.
 *
 *  The complete API of the Poisson distribution is described in poisson.h.
 *  The device API expects Poisson distribution objects to be stored into
 *  global memory by default.
 *  See @ref mem_types for information on using other types of memory.
 *
 *  **WARNING:** The complementary distribution function, which can be
 *  evaluated with clprobdistPoissonComplCDF() and
 *  clprobdistPoissonComplCDFWithObject(),
 *  is defined as \f$\bar F(j) = \mathbb P[X \geq j]\f$ (for integers \f$j\f$),
 *  so that for the Poisson distribution clProbDist,
 *  \f$F(j) + \bar F(j) \neq 1\f$ since both include the term
 *  \f$\mathbb P[X = j]\f$.
 */

/*! @example DocsTutorial/example1.c
 *
 *  Complete code for @ref example_ci
 */

/*! @example DocsTutorial/example2.c
 *
 *  Complete code for @ref example_poissongen
 */

/*! @example DocsTutorial/example3.c
 *
 *  Complete code for @ref example_poissongen_obj
 */

/*! @example DocsTutorial/example4.c
 *
 *  Host code for @ref example_poissongen_device
 */

/*! @example DocsTutorial/example4_kernel.cl
 *
 *  Device code for @ref example_poissongen_device
 */


/*! @file clProbDist.h
*  @brief Library definitions common to all probability distributions
*  @see clProbDist_template.h
*/

#pragma once
#ifndef CLPROBDIST_H
#define CLPROBDIST_H

#ifdef CLPROBDIST_SINGLE_PRECISION
#error "CLPROBDIST_SINGLE_PRECISION option not yet implemented"
#endif

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#if defined ( WIN32 )
#define __func__ __FUNCTION__
#endif

#define constant const

/*! @brief Error codes
*
*  Most library functions return an error status indicating the success or
*  error state of the operation carried by the function.
*  In case of success, the error status is set to `CL_SUCCESS`.
*  Otherwise, an error message can be retrieved by invoking
*  clprobdistGetErrorString().
*
*  @note In naming this type clprobdistStatus, we follow the convention from clFFT
*  and clBLAS, where the homologous types are name clfftStatus and
*  clblasStatus, respectively.
*/
typedef enum clprobdistStatus_ {
	CLPROBDIST_SUCCESS = CL_SUCCESS,
	CLPROBDIST_OUT_OF_RESOURCES = CL_OUT_OF_RESOURCES,
	CLPROBDIST_INVALID_VALUE = CL_INVALID_VALUE,
	/* ... */
	CLPROBDIST_INVALID_ENVIRONMENT,
	CLPROBDIST_FUNCTION_NOT_IMPLEMENTED,
	CLPROBDIST_ABSTRACT_FUNCTION
} clprobdistStatus;

#ifdef __cplusplus
extern "C" {
#endif

	/*! @brief Retrieve the last error message.
	*
	*  The buffer containing the error message is internally allocated and must
	*  not be freed by the client.
	*
	*  @return     Error message or `NULL`.
	*/
	const char* clprobdistGetErrorString();

	/*! @brief Generate an include option string for use with the OpenCL C compiler
	*
	*  Generate and return "-I${CLPROBDIST_ROOT}/include", where \c ${CLPROBDIST_ROOT} is
	*  the value of the \c CLPROBDIST_ROOT environment variable.
	*  This string is meant to be passed as an option to the OpenCL C compiler for
	*  programs that make use of the clProbDist device-side headers.
	*  If the \c CLPROBDIST_ROOT environment variable is not defined, it
        *  defaults to `/usr` if the file
        *  `/usr/include/clProbDist/clProbDist.h` exists, else to the current
	*  directory of execution of the program.
	*
	*  A static buffer is return and need not be released; it could change upon
	*  successive calls to the function.
	*
	*  An error is returned in \c err if the preallocated buffer is too small to
	*  contain the include string.
	*
	*  @param[out]     err         Error status variable, or `NULL`.
	*
	*  @return An OpenCL C compiler option to indicate where to find the
	*  device-side clProbDist headers.
	*/
	const char* clprobdistGetLibraryDeviceIncludes(cl_int* err);

	/*! @brief Retrieve the library installation path
	*
        *  @return Value of the CLPROBDIST_ROOT environment variable, if
        *  defined; else, `/usr` if the file
        *  `/usr/include/clProbDist/clProbDist.h` exists; or, the current
	*  directory (.) of execution of the program otherwise.
	*/
	const char* clprobdistGetLibraryRoot();

#ifdef __cplusplus
}
#endif

#endif /* CLPROBDIST_H *
* vim: syntax=c.doxygen spell spelllang=en fdm=syntax fdls=0 expandtab
*/
