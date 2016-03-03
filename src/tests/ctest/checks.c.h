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
 *
 ***********************************************************************
 Copyright (c) 2015 Advanced Micro Devices, Inc. 
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without 
 modification, are permitted provided that the following conditions 
 are met:
 
 1. Redistributions of source code must retain the above copyright 
 notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright 
 notice, this list of conditions and the following disclaimer in the 
 documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 ***********************************************************************
 */

/*! @file checks.c.h
 *  @brief Tests that do not depend on the floating-point precision.
 *
 *  These tests must be compiled once for every generator.
 */

#include "mangle.h"
#include "util.h"
#include "../client/common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifndef CTEST_REL_TOL
#define CTEST_REL_TOL 5.0 * DBL_EPSILON
#endif
#define SUCCESS_STR(ret) ((ret) == EXIT_SUCCESS ? "  SUCCESS" : "* FAILURE")
#define DBL_EQUALS(a, b) (fabs((a) - (b)) < (CTEST_REL_TOL) * fmax(fmax(fabs(a), fabs(b)), 1.0))

#ifdef CTEST_PROBDIST_DISCRETE
typedef cl_int var_type;
#define clprobdistProbOrDensity             clprobdistProb
#define clprobdistProbOrDensityWithObject   clprobdistProbWithObject
#else
typedef cl_double var_type;
#define clprobdistProbOrDensity             clprobdistDensity
#define clprobdistProbOrDensityWithObject   clprobdistDensityWithObject
#endif


/*! @brief Compares density, CDF, InverseCDF and ComplCDF values with expected
 * values on the host.
 */
int CTEST_MANGLE(checkHost)()
{
  size_t num_params = CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams));
  size_t num_values = CTEST_ARRAY_SIZE(CTEST_MANGLE(inputCDFValues));
  clprobdistStatus err;
  size_t test_count = 0;
  size_t failed_count = 0;

  for (size_t k = 0; k < num_params; k++) {

    clprobdistObject* dist = clprobdistCreate(
        CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]),
        NULL, &err);
    check_error(err, NULL);

    for (size_t i = 0; i < num_values; i++) {

      cl_double u                    = CTEST_MANGLE(inputCDFValues)[i];
      cl_double expectedInverseCDF   = CTEST_MANGLE(expectedInverseCDFValues)[k][i];
      cl_double observedInverseCDF   = clprobdistInverseCDF(CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]), u, &err); check_error(err, NULL);
      cl_double observedInverseCDF_o = clprobdistInverseCDFWithObject(dist, u, &err); check_error(err, NULL);
      var_type  x                    = observedInverseCDF_o;
      cl_double expectedDensity      = CTEST_MANGLE(expectedDensityValues)[k][i];
      cl_double observedDensity      = clprobdistProbOrDensity(CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]), x, &err); check_error(err, NULL);
      cl_double observedDensity_o    = clprobdistProbOrDensityWithObject(dist, x, &err); check_error(err, NULL);
#ifdef CTEST_PROBDIST_DISCRETE
      cl_double expectedCDF          = CTEST_MANGLE(expectedCDFValues)[k][i];
#else
      cl_double expectedCDF          = u;
#endif
      cl_double observedCDF          = clprobdistCDF(CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]), x, &err); check_error(err, NULL);
      cl_double observedCDF_o        = clprobdistCDFWithObject(dist, x, &err); check_error(err, NULL);
      cl_double expectedComplCDF     = 1.0 - expectedCDF;
#ifdef CTEST_PROBDIST_DISCRETE
      expectedComplCDF               += expectedDensity;
#endif
      cl_double observedComplCDF     = clprobdistComplCDF(CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]), x, &err); check_error(err, NULL);
      cl_double observedComplCDF_o   = clprobdistComplCDFWithObject(dist, x, &err); check_error(err, NULL);

      test_count++;

      if (
          !DBL_EQUALS(observedInverseCDF,    expectedInverseCDF)  ||
          !DBL_EQUALS(observedInverseCDF_o,  expectedInverseCDF)  ||
          !DBL_EQUALS(observedDensity,       expectedDensity)     ||
          !DBL_EQUALS(observedDensity_o,     expectedDensity)     ||
          !DBL_EQUALS(observedCDF,           expectedCDF)         ||
          !DBL_EQUALS(observedCDF_o,         expectedCDF)         ||
          !DBL_EQUALS(observedComplCDF,      expectedComplCDF)    ||
          !DBL_EQUALS(observedComplCDF_o,    expectedComplCDF)
         ) {

        failed_count++;

        if (ctestVerbose) {

          printf("\n%4sValues do not match at params", "");

#if CTEST_PROBDIST_PARAM_COUNT >= 1
            printf(" %g", CTEST_MANGLE(inputDistParams)[k].p1);
#endif
#if CTEST_PROBDIST_PARAM_COUNT >= 2
            printf(" %g", CTEST_MANGLE(inputDistParams)[k].p2);
#endif
#if CTEST_PROBDIST_PARAM_COUNT >= 3
            printf(" %d", CTEST_MANGLE(inputDistParams)[k].p3);
#endif
          printf(" and u = %g.\n", CTEST_MANGLE(inputCDFValues)[i]);

          if (!DBL_EQUALS(observedInverseCDF, expectedInverseCDF) ||
              !DBL_EQUALS(observedInverseCDF_o, expectedInverseCDF)) {
            printf("\n%8sclprobdistInverseCDF() / clprobdistInverseCDFWithObject() / expected:\n", "");
            printf("%12s%.18g / %.18g / %.18g\n", "", observedInverseCDF, observedInverseCDF_o, expectedInverseCDF);
          }
          if (!DBL_EQUALS(observedDensity, expectedDensity) ||
              !DBL_EQUALS(observedDensity_o, expectedDensity)) {
            printf("\n%8sclprobdistDensity() / clprobdistDensityWithObject() / expected:\n", "");
            printf("%12s%.18g / %.18g / %.18g\n", "", observedDensity, observedDensity_o, expectedDensity);
          }
          if (!DBL_EQUALS(observedCDF, expectedCDF) ||
              !DBL_EQUALS(observedCDF_o, expectedCDF)) {
            printf("\n%8sclprobdistCDF() / clprobdistCDFWithObject() / expected:\n", "");
            printf("%12s%.18g / %.18g / %.18g\n", "", observedCDF, observedCDF_o, expectedCDF);
          }
          if (!DBL_EQUALS(observedComplCDF, expectedComplCDF) ||
              !DBL_EQUALS(observedComplCDF_o, expectedComplCDF)) {
            printf("\n%8sclprobdistComplCDF() / clprobdistComplCDFWithObject() / expected:\n", "");
            printf("%12s%.18g / %.18g / %.18g\n", "", observedComplCDF, observedComplCDF_o, expectedComplCDF);
          }
          printf("\n");
        }
      }
    }

    err = clprobdistDestroy(dist);
    check_error(err, NULL);
  }

  int ret = failed_count == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

  char subreport[128];
  if (failed_count > 0)
    sprintf(subreport, "%6" SIZE_T_FORMAT " / %6" SIZE_T_FORMAT, failed_count, test_count);
  else
    sprintf(subreport, "%15" SIZE_T_FORMAT, test_count);

  printf("%s  %s tests  -  %16s Host API\n",
      SUCCESS_STR(ret), subreport, PROBDIST_TYPE_S);

  return ret;
}


#if 0
int CTEST_MANGLE(plot)()
{
  int ret = EXIT_SUCCESS;
  size_t num_params = CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams));
  clprobdistStatus err;

  for (size_t k = 0; k < num_params && ret == EXIT_SUCCESS; k++) {

    printf("\nParams");
    for (size_t j = 0; j < CTEST_PROBDIST_PARAM_COUNT; j++)
      printf(" %g", CTEST_MANGLE(inputDistParams)[k][j]);
    printf(":\n");

    printf("%12s\t%12s\t%12s\t%12s\n", "x", "density", "cdf", "ccdf");

    clprobdistObject* dist = clprobdistCreate(
        CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]),
        NULL, &err);
    check_error(err, NULL);

    size_t n = 30;
    for (size_t i = 0; i < n; i++) {

      cl_double u       = (i + 0.5) / n;
      cl_double x       = clprobdistInverseCDFWithObject(dist, u, &err);
      cl_double density = clprobdistDensityWithObject(dist, x, &err); check_error(err, NULL);
      cl_double cdf     = clprobdistCDFWithObject(dist, x, &err); check_error(err, NULL);
      cl_double ccdf    = clprobdistComplCDFWithObject(dist, x, &err); check_error(err, NULL);

      printf("%12g\t%12g\t%12f\t%12f\n", x, density, cdf, ccdf);
    }

    err = clprobdistDestroy(dist);
    check_error(err, NULL);
  }

  return ret;
}
#endif

/*! @brief Helper function for checkDeviceOperations()
 */
static int hostOperations(
    size_t            gsize,
    size_t            quota,
    clprobdistObject* dist,
    double*           prob,
    double*           cdf,
    double*           compl_cdf,
    var_type*         inverse_cdf)
{
  for (size_t i = 0; i < quota; i++) {
    for (size_t gid = 0; gid < gsize; gid++) {
      cl_uint k      = i * gsize + gid;
      cl_double u    = ((double) k) / (gsize * quota);
      var_type x     = clprobdistInverseCDFWithObject(dist, u, NULL);
#ifdef CTEST_PROBDIST_DISCRETE
      prob[k]        = clprobdistProbWithObject(dist, x, NULL);
#else
      prob[k]        = clprobdistDensityWithObject(dist, x, NULL);
#endif
      cdf[k]         = clprobdistCDFWithObject(dist, x, NULL);
      compl_cdf[k]   = clprobdistComplCDFWithObject(dist, x, NULL);
      inverse_cdf[k] = x;
    }
  }
  return EXIT_SUCCESS;
}

/*! @brief Structure for use with checkDevice()
 */
typedef struct DeviceParams_ {
  size_t            num_work_items;
  cl_uint           quota_per_work_item;
  DistParams*       dist_params;
  size_t            dist_size;
  clprobdistObject* dist;
  MemType           mem_type;  // defined in mangle.h
  double*           out_prob;
  double*           out_cdf;
  double*           out_compl_cdf;
  var_type*         out_inverse_cdf;
} DeviceParams;

/*! @brief Helper function for checkDevice()
*/
static int deviceOperations(cl_context context, cl_device_id device, cl_command_queue queue, void* data_)
{
  cl_int err;
  const DeviceParams* data = (const DeviceParams*) data_;
  size_t output_count = data->num_work_items * data->quota_per_work_item;

  cl_mem dist_buf;
  if (data->mem_type == MEM_PARAMS) {
    dist_buf = clCreateBuffer(context,
      CL_MEM_HOST_WRITE_ONLY | CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
      sizeof(DistParams), data->dist_params, &err);
  }
  else {
    dist_buf = clCreateBuffer(context,
      CL_MEM_HOST_WRITE_ONLY | CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
      data->dist_size, data->dist, &err);
  }
  check_error(err, "cannot create dist buffer");
  cl_mem prob_buf = clCreateBuffer(context,
      CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, 
      output_count * sizeof(double), NULL, &err);
  check_error(err, "cannot create prob buffer");
  cl_mem cdf_buf = clCreateBuffer(context,
      CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, 
      output_count * sizeof(double), NULL, &err);
  check_error(err, "cannot create cdf buffer");
  cl_mem compl_cdf_buf = clCreateBuffer(context,
      CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, 
      output_count * sizeof(double), NULL, &err);
  check_error(err, "cannot create compl_cdf buffer");
  cl_mem inverse_cdf_buf = clCreateBuffer(context,
      CL_MEM_HOST_READ_ONLY | CL_MEM_WRITE_ONLY, 
      output_count * sizeof(var_type), NULL, &err);
  check_error(err, "cannot create inverse_cdf buffer");

  const char* func_suffix = data->mem_type == MEM_PARAMS ? "" : "WithObject";
  char params_arg1[256] = "\000";
  char params_arg2[256] = "\000";
  char params_str[256]  = "\000";
  char sync_block[512]  = "\000";
  if (data->mem_type == MEM_PARAMS) {
    for (int j = 0; j < CTEST_PROBDIST_PARAM_COUNT; j++) {
      strcpy(params_arg1, "                   const __global DistParams* params,\n");
      char buf[100];
      sprintf(buf, "%sparams->p%d", j > 0 ? ", " : "", j + 1);
      strcat(params_str, buf);
    }
  }
  else {
    strcpy(params_str, "dist");
    sprintf(params_arg1,
        "                   %s const clprobdist" PROBDIST_TYPE_S "* %sdist, \n",
        data->mem_type == MEM_CONSTANT ? "__constant" : "__global",
        data->mem_type == MEM_LOCAL || data->mem_type == MEM_PRIVATE ? "g" : "");
    if (data->mem_type == MEM_LOCAL) {
      strcpy(params_arg2,
          "                   __local  clprobdist" PROBDIST_TYPE_S "* dist, \n");
      strcpy(sync_block,
          "    int group_size = get_local_size(0);\n"
          "    if (gid % group_size == 0)\n"
          "        clprobdist" PROBDIST_TYPE_S "CopyOverFromGlobal(dist, gdist);\n"
          "    barrier(CLK_LOCAL_MEM_FENCE);\n");
    }
    else if (data->mem_type == MEM_PRIVATE) {
      strcpy(params_str, "&pdist");
      strcpy(sync_block,
          "    clprobdist" PROBDIST_TYPE_S " pdist;\n"
          "    clprobdist" PROBDIST_TYPE_S "CopyOverFromGlobal(&pdist, gdist);\n");
    }
  }

  char source[1500];
  sprintf(source,
      //"#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n"
      "%s\n"
      "#include " PROBDIST_DEVICE_HEADER_S "\n"
#ifdef CTEST_PROBDIST_DISCRETE
      "typedef int    var_type;\n"
#else
      "typedef double var_type;\n"
#endif
      "typedef struct { " _CTEST_STR_(CTEST_PROBDIST_PARAM_DEF) " } DistParams;\n"
      "__kernel void test(\n"
      "                   uint quota,\n"
      "%s"
      "%s"
      "                   __global double* prob,\n"
      "                   __global double* cdf, \n"
      "                   __global double* compl_cdf, \n"
      "                   __global var_type* inverse_cdf) {\n"
      "    int gid = get_global_id(0);\n"
      "    int gsize = get_global_size(0);\n"
      "%s"
      "    for (uint i = 0; i < quota; i++) {\n"
      "        uint k         = i * gsize + gid;\n"
      "        double u       = ((double) k) / (gsize * quota);\n"
      "        var_type x     = clprobdist" PROBDIST_TYPE_S "InverseCDF%s(%s, u, (void*)0);\n"
#ifdef CTEST_PROBDIST_DISCRETE
      "        prob[k]        = clprobdist" PROBDIST_TYPE_S "Prob%s(%s, x, (void*)0);\n"
#else
      "        prob[k]        = clprobdist" PROBDIST_TYPE_S "Density%s(%s, x, (void*)0);\n"
#endif
      "        cdf[k]         = clprobdist" PROBDIST_TYPE_S "CDF%s(%s, x, (void*)0);\n"
      "        compl_cdf[k]   = clprobdist" PROBDIST_TYPE_S "ComplCDF%s(%s, x, (void*)0);\n"
      "        inverse_cdf[k] = x;\n"
      "    }\n"
      "}\n",
      data->mem_type == MEM_PRIVATE  ? "#define CLPROBDIST_" PROBDIST_TYPE_SC "_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE"  :
      data->mem_type == MEM_LOCAL    ? "#define CLPROBDIST_" PROBDIST_TYPE_SC "_OBJ_MEM CLPROBDIST_MEM_TYPE_LOCAL"    :
      data->mem_type == MEM_CONSTANT ? "#define CLPROBDIST_" PROBDIST_TYPE_SC "_OBJ_MEM CLPROBDIST_MEM_TYPE_CONSTANT" :
      data->mem_type == MEM_GLOBAL   ? "#define CLPROBDIST_" PROBDIST_TYPE_SC "_OBJ_MEM CLPROBDIST_MEM_TYPE_GLOBAL"   : "",
    params_arg1, params_arg2, sync_block,
    func_suffix, params_str,
    func_suffix, params_str,
    func_suffix, params_str,
    func_suffix, params_str);

  if (ctestVerbose >= 2)
    printf("source:\n\n%s\nsource length:%" SIZE_T_FORMAT "\n", source, strlen(source));

  const char* sources[] = { source };
  cl_program program = clCreateProgramWithSource(context, 1, sources, NULL, &err);
  check_error(err, "cannot create program");

  const char* includes = clprobdistGetLibraryDeviceIncludes(&err);
  check_error(err, NULL);

  err = clBuildProgram(program, 0, NULL, includes, NULL, NULL);
  if (err < 0)
    write_build_log(stderr, program, device);
  check_error(err, "cannot build program");

  cl_kernel kernel = clCreateKernel(program, "test", &err);
  check_error(err, "cannot create kernel");

  size_t iarg = 0;
  err  = clSetKernelArg(kernel, iarg++, sizeof(data->quota_per_work_item),  &data->quota_per_work_item);
  err |= clSetKernelArg(kernel, iarg++, sizeof(dist_buf),                   &dist_buf);
  if (data->mem_type == MEM_LOCAL)
    // create local memory buffer for distribution object
    err |= clSetKernelArg(kernel, iarg++, data->dist_size,                  NULL);
  err |= clSetKernelArg(kernel, iarg++, sizeof(prob_buf),                   &prob_buf);
  err |= clSetKernelArg(kernel, iarg++, sizeof(cdf_buf),                    &cdf_buf);
  err |= clSetKernelArg(kernel, iarg++, sizeof(compl_cdf_buf),              &compl_cdf_buf);
  err |= clSetKernelArg(kernel, iarg++, sizeof(inverse_cdf_buf),            &inverse_cdf_buf);
  check_error(err, "cannot set kernel arguments");

  cl_event ev;
  err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &data->num_work_items, NULL, 0, NULL, &ev);
  check_error(err, "cannot enqueue kernel");

  err = clWaitForEvents(1, &ev);
  check_error(err, "error waiting for events");

  // read output values
  err = clEnqueueReadBuffer(queue, prob_buf,        CL_TRUE, 0, output_count * sizeof(double),   data->out_prob,        0, NULL, NULL);
  check_error(err, "cannot read prob buffer");
  err = clEnqueueReadBuffer(queue, cdf_buf,         CL_TRUE, 0, output_count * sizeof(double),   data->out_cdf,         0, NULL, NULL);
  check_error(err, "cannot read cdf buffer");
  err = clEnqueueReadBuffer(queue, compl_cdf_buf,   CL_TRUE, 0, output_count * sizeof(double),   data->out_compl_cdf,   0, NULL, NULL);
  check_error(err, "cannot read compl_cdf buffer");
  err = clEnqueueReadBuffer(queue, inverse_cdf_buf, CL_TRUE, 0, output_count * sizeof(var_type), data->out_inverse_cdf, 0, NULL, NULL);
  check_error(err, "cannot read inverse_cdf buffer");

  clReleaseEvent(ev);
  clReleaseMemObject(dist_buf);
  clReleaseMemObject(prob_buf);
  clReleaseMemObject(cdf_buf);
  clReleaseMemObject(compl_cdf_buf);
  clReleaseMemObject(inverse_cdf_buf);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  return EXIT_SUCCESS;
}

/*! @brief Device operations.
 */
int CTEST_MANGLE(checkDeviceHelper)(
  size_t  num_work_items,
  cl_uint quota_per_work_item,
  MemType mem_type,
  const DeviceSelect* dev)
{
  size_t output_count = num_work_items * quota_per_work_item;

  double*   device_prob        = (double*)   malloc(output_count * sizeof(double));
  double*   device_cdf         = (double*)   malloc(output_count * sizeof(double));
  double*   device_compl_cdf   = (double*)   malloc(output_count * sizeof(double));
  var_type* device_inverse_cdf = (var_type*) malloc(output_count * sizeof(double));
  double*   host_prob          = (double*)   malloc(output_count * sizeof(double));
  double*   host_cdf           = (double*)   malloc(output_count * sizeof(double));
  double*   host_compl_cdf     = (double*)   malloc(output_count * sizeof(double));
  var_type* host_inverse_cdf   = (var_type*) malloc(output_count * sizeof(double));

  size_t num_params = CTEST_ARRAY_SIZE(CTEST_MANGLE(inputDistParams));

  clprobdistStatus err;
  size_t test_count = 0;
  size_t failed_count = 0;

  for (size_t k = 0; k < num_params; k++) {

    size_t dist_size;
    clprobdistObject* dist = clprobdistCreate(
        CTEST_DIST_PARAMS(CTEST_MANGLE(inputDistParams)[k]),
        &dist_size,
        &err);
    check_error(err, NULL);

    hostOperations(
        num_work_items,
        quota_per_work_item,
        dist,
        host_prob,
        host_cdf,
        host_compl_cdf,
        host_inverse_cdf);

    DeviceParams params = {
      num_work_items,
      quota_per_work_item,
      &CTEST_MANGLE(inputDistParams)[k],
      dist_size,
      dist,
      mem_type,
      device_prob,
      device_cdf,
      device_compl_cdf,
      device_inverse_cdf
    };

    call_with_opencl(dev->platform_index, dev->device_type, dev->device_index, &deviceOperations, &params, false);

    test_count += 4 * output_count;

    // compare output values
    for (size_t i = 0; i < output_count; i++) {

      if (!DBL_EQUALS(device_prob[i], host_prob[i])) {
        failed_count++;
        if (ctestVerbose) {
          printf("\n%8sdevice / host prob [%" SIZE_T_FORMAT "]:\n", "", i);
          printf("%12s%.18g / %.18g\n", "", device_prob[i], host_prob[i]);
        }
      }
      if (!DBL_EQUALS(device_cdf[i], host_cdf[i])) {
        failed_count++;
        if (ctestVerbose) {
          printf("\n%8sdevice / host cdf [%" SIZE_T_FORMAT "]:\n", "", i);
          printf("%12s%.18g / %.18g\n", "", device_cdf[i], host_cdf[i]);
        }
      }
      if (!DBL_EQUALS(device_compl_cdf[i], host_compl_cdf[i])) {
        failed_count++;
        if (ctestVerbose) {
          printf("\n%8sdevice / host compl_cdf [%" SIZE_T_FORMAT "]:\n", "", i);
          printf("%12s%.18g / %.18g\n", "", device_compl_cdf[i], host_compl_cdf[i]);
        }
      }
      if (!DBL_EQUALS(device_inverse_cdf[i], host_inverse_cdf[i])) {
        failed_count++;
        if (ctestVerbose) {
          printf("\n%8sdevice / host inverse_cdf [%" SIZE_T_FORMAT "]:\n", "", i);
          printf("%12s%.18g / %.18g\n", "", (cl_double) device_inverse_cdf[i], (cl_double) host_inverse_cdf[i]);
        }
      }
    }

    err = clprobdistDestroy(dist);
    check_error(err, NULL);
  }

#if 0
  if (ctestVerbose >= 2) {
    printf("prob:\n");
    write_array(stdout, 4, "%16.12f", output_count, device_prob);
    printf("cdf:\n");
    write_array(stdout, 4, "%16.12f", output_count, device_cdf);
    printf("compl_cdf:\n");
    write_array(stdout, 4, "%16.12f", output_count, device_compl_cdf);
    printf("inverse_cdf:\n");
    write_array(stdout, 4, "%16.12f", output_count, device_inverse_cdf);
  }
#endif

  free(device_prob);
  free(device_cdf);
  free(device_compl_cdf);
  free(device_inverse_cdf);
  free(host_prob);
  free(host_cdf);
  free(host_compl_cdf);
  free(host_inverse_cdf);

  int ret = failed_count == 0 ? EXIT_SUCCESS : EXIT_FAILURE;

  char subreport[128];
  if (failed_count > 0)
    sprintf(subreport, "%6" SIZE_T_FORMAT " / %6" SIZE_T_FORMAT, failed_count, test_count);
  else
    sprintf(subreport, "%15" SIZE_T_FORMAT, test_count);

  printf("%s  %s tests  -  %16s Device API with %s (%" SIZE_T_FORMAT ",%d)\n",
      SUCCESS_STR(ret), subreport, PROBDIST_TYPE_S,
      mem_type == MEM_PRIVATE  ? "private objects"  :
      mem_type == MEM_LOCAL    ? "local objects"    :
      mem_type == MEM_CONSTANT ? "constant objects" :
      mem_type == MEM_GLOBAL   ? "global objects"   : "parameters",
      num_work_items, quota_per_work_item);


  return ret;
}

int CTEST_MANGLE(checkDevice)(const DeviceSelect* dev)
{
  size_t  num_work_items        = 1 << 10;
  cl_uint quota_per_work_item   = 1 << 5;

  int ret = EXIT_SUCCESS;

  size_t num_mem_types = CTEST_ARRAY_SIZE(CTEST_MANGLE(objMemTypes));
  for (size_t i = 0; i < num_mem_types; i++)
    ret |= CTEST_MANGLE(checkDeviceHelper)(
        num_work_items,
        quota_per_work_item,
        CTEST_MANGLE(objMemTypes)[i],
        dev);

  return ret;
}
