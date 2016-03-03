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
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 */

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include "../common.h"

#include <clProbDist/poisson.h>
#include <clRNG/mrg31k3p.h>

int task(cl_context context, cl_device_id device, cl_command_queue queue, void* data)
{
    cl_int err;
    size_t numWorkItems = *(size_t*)data;

    /*! [create distribution object] */
    size_t dist_buf_size;
    clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, &dist_buf_size, NULL);
    /*! [create distribution object] */
    check_error(err, NULL);
    /*! [create streams] */
    size_t streams_buf_size;
    clrngMrg31k3pStream* streams = clrngMrg31k3pCreateStreams(NULL, numWorkItems,
                               &streams_buf_size, (clrngStatus *)&err);
    /*! [create streams] */
    check_error(err, "cannot create random streams array");

    cl_program program = build_program_from_file(context, device,
	    "client/DocsTutorial/example4_kernel.cl", NULL, PATH_RELATIVE_TO_LIB,
	    clrngGetLibraryDeviceIncludes(&err));
    check_error(err, NULL);
    cl_kernel kernel = clCreateKernel(program, "my_kernel", &err);
    check_error(err, "cannot create kernel");


    // Create buffer to transfer streams to the device.
    /*! [create streams buffer] */
    cl_mem streams_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
                                    streams_buf_size, streams, &err);
    /*! [create streams buffer] */
    check_error(err, "cannot create streams buffer");
    /*! [create dist buffer] */
    // Create buffer to transfer streams to the device.
    cl_mem dist_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
                                    dist_buf_size, dist, &err);
    /*! [create dist buffer] */
    check_error(err, "cannot create dist buffer");
    // Create buffer to transfer output back from the device.
    /*! [create output buffer] */
    cl_mem out_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, 
                                     numWorkItems * sizeof(cl_double), NULL, &err);
    /*! [create output buffer] */
    check_error(err, "cannot create output buffer");

    // Set kernel arguments
    err  = clSetKernelArg(kernel, 0, sizeof(dist_buf),    &dist_buf);
    err |= clSetKernelArg(kernel, 1, sizeof(streams_buf), &streams_buf);
    err |= clSetKernelArg(kernel, 2, sizeof(out_buf),     &out_buf);
    check_error(err, "cannot set kernel arguments");

    // Enqueue the kernel on device.
    cl_event ev;
    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &numWorkItems, NULL, 0, NULL, &ev);
    check_error(err, "cannot enqueue kernel");

    // Wait for all work items to finish.
    err = clWaitForEvents(1, &ev);
    check_error(err, "error waiting for events");

    cl_ulong t0, t1;
    clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(t0), &t0, NULL);
    clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END,   sizeof(t1), &t1, NULL);
    cl_ulong total_time = t1 - t0;

    // Retrieve the contents of the output buffer from the device.
    cl_double* out = (cl_double*) malloc(numWorkItems * sizeof(cl_double));
    err = clEnqueueReadBuffer(queue, out_buf, CL_TRUE, 0,
	    numWorkItems * sizeof(out[0]), out, 0, NULL, NULL);
    check_error(err, "cannot read output buffer");

    // printf("output buffer:\n");
    // for (int j = 0; j < numWorkItems; j++)
    //     printf("    %f\n", out[j]);

    double out_sum = 0.0;
    for (int j = 0; j < numWorkItems; j++)
	out_sum += out[j];
    printf("\naverage output: %.3f\n", out_sum / numWorkItems);

    printf("\nprocessing time: %1.2f\n", total_time * 1e-9);

    clprobdistPoissonDestroy(dist);
    clrngMrg31k3pDestroyStreams(streams);
    free(out);
    clReleaseEvent(ev);
    clReleaseMemObject(dist_buf);
    clReleaseMemObject(streams_buf);
    clReleaseMemObject(out_buf);
    clReleaseKernel(kernel);
    clReleaseProgram(program);

    return EXIT_SUCCESS;
}


int main()
{
    size_t numWorkItems = 1024;
    return call_with_opencl(0, CL_DEVICE_TYPE_CPU, 0, &task, &numWorkItems, true);
}


