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
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <clProbDist/clProbDist.h>
#include <clProbDist/exponential.h>


int main( void )
{
    cl_int err;
    cl_platform_id platform = 0;
    cl_device_id device = 0;
    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
    cl_context ctx = 0;
    cl_command_queue queue = 0;
    cl_program program = 0;
    cl_kernel kernel = 0;
    cl_event event = 0;
    cl_mem bufIn, bufOut;
    double *out;
    const char *includes;
    char buildLog[4096];
    size_t numWorkItems = 64;
    size_t i;
    clprobdistExponential *dist = 0;
    size_t distBufferSize = 0;
    size_t kernelLines = 0;

    /* Sample kernel that calls clProbDist device-side interfaces to compute the
       inverse CDF of the exponential distribution at different quantiles */
    const char *kernelSrc[] = {
        "#include <clProbDist/exponential.clh>                             \n",
        "                                                                  \n",
        "__kernel void example(__global const clprobdistExponential *dist, \n",
        "                      __global double *out)                       \n",
        "{                                                                 \n",
        "    int gid = get_global_id(0);                                   \n",
        "    int gsize = get_global_size(0);                               \n",
        "    double quantile = (gid + 0.5) / gsize;                        \n",
        "    out[gid] = clprobdistExponentialInverseCDFWithObject(         \n",
        "                                      dist, quantile, (void *)0); \n",
        "}                                                                 \n",
    };

    /* Setup OpenCL environment. */
    err = clGetPlatformIDs( 1, &platform, NULL );
    err = clGetDeviceIDs( platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL ); // FIXME
    err = clGetDeviceIDs( platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL );

    props[1] = (cl_context_properties)platform;
    ctx = clCreateContext( props, 1, &device, NULL, NULL, &err );
    #ifdef CL_VERSION_2_0
        queue = clCreateCommandQueueWithProperties( ctx, device, (cl_queue_properties[]){0}, &err );
    #else
        queue = clCreateCommandQueue( ctx, device, 0, &err );
    #endif


    /* Make sure CLPROBDIST_ROOT is specified to get library path */
    includes = clprobdistGetLibraryDeviceIncludes(&err);
    if(err != CL_SUCCESS) printf("\n%s\nSpecify environment variable CLPROBDIST_ROOT as described\n", clprobdistGetErrorString());

    /* Create sample kernel */
    kernelLines = sizeof(kernelSrc) / sizeof(kernelSrc[0]);
    program = clCreateProgramWithSource(ctx, kernelLines, kernelSrc, NULL, &err);
    err = clBuildProgram(program, 1, &device, includes, NULL, NULL);
    if(err != CL_SUCCESS)
    {
        printf("\nclBuildProgram has failed\n");
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 4096, buildLog, NULL);
        printf("%s", buildLog);
        exit(1);
    }
    kernel = clCreateKernel(program, "example", &err);

    /* Create an exponential distribution with unit mean */
    dist = clprobdistExponentialCreate(1.0, &distBufferSize, (clprobdistStatus *)&err);

    /* Create buffers for the kernel */
    bufIn = clCreateBuffer(ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, distBufferSize, dist, &err);
    bufOut = clCreateBuffer(ctx, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, numWorkItems * sizeof(cl_double), NULL, &err);

    /* Setup the kernel */
    err = clSetKernelArg(kernel, 0, sizeof(bufIn),  &bufIn);
    err = clSetKernelArg(kernel, 1, sizeof(bufOut), &bufOut);

    /* Execute the kernel and read back results */
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &numWorkItems, NULL, 0, NULL, &event);
    err = clWaitForEvents(1, &event);
    out = (double *)malloc(numWorkItems * sizeof(out[0]));
    err = clEnqueueReadBuffer(queue, bufOut, CL_TRUE, 0, numWorkItems * sizeof(out[0]), out, 0, NULL, NULL);

    /* Display results and check that the original quantile is recovered by evaluating the CDF */
    for (i = 0; i < numWorkItems; i++)
        printf("quantile %.3f :   %9.3e   ->   %.3f\n", (i + 0.5) / numWorkItems, out[i], clprobdistExponentialCDFWithObject(dist, out[i], &err));

    /* Release allocated resources */
    clReleaseEvent(event);
    free(out);
    clReleaseMemObject(bufIn);
    clReleaseMemObject(bufOut);

    clReleaseKernel(kernel);
    clReleaseProgram(program);

    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);

    return 0;
}
