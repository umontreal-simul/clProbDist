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
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "Options.h"

#include "../common.h"
#include <clRNG/mrg32k3a.h>

#include <clProbDist/clProbDist.h>
#include <clProbDist/poisson.h>
#include <clProbDist/normal.h>

static char _options[255];
char* BuildOption(char _options[255], int n2, double lambda)
{
	char _n2[7], _lambda[20];

	sprintf(_n2, "%d", n2);
	sprintf(_lambda, "%f", lambda);

	cl_int err;
	strcpy(_options, clrngGetLibraryDeviceIncludes(&err));
	check_error(err, "cannot get clRNG device-side headers");

	strcat(_options, " -Dparam_lambda=\"");
	strcat(_options, _lambda);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_n2=\"");
	strcat(_options, _n2);

	return strcat(_options, "\"");
}

void performanceOnCPU_Poisson(int n, double lambda, clrngMrg32k3aStream* stream, ExecType execType)
{
	//Create Poisson Dist
	clprobdistStatus err;
	size_t poissonBufSize = 0;
	clprobdistPoisson* poissonDist = NULL;

	poissonDist = clprobdistPoissonCreate(lambda, &poissonBufSize, &err);
	check_error(err, "%s(): cannot create Poisson distribution", __func__);

	if (execType == simWithObject){
		//*******************************************************
		//call the Poisson Object n times
		clock_t start = clock();
		for (int i = 0; i < n; i++){

			double u = clrngMrg32k3aRandomU01(stream);
			clprobdistPoissonInverseCDFWithObject(poissonDist, u, NULL);

		}

		//Compute Execution Time
		clock_t end = clock();
		float CPU_time = (float)(end - start) / CLOCKS_PER_SEC;
		printf("\n\tLambda : %1.1f (O) : Total CPU time (sec.) : %1.6f\n", lambda, CPU_time);
	}
	else
	{
		//*******************************************************
		//call the Poisson static function n times
		clock_t start = clock();
		for (int i = 0; i < n; i++) {
			double u = clrngMrg32k3aRandomU01(stream);
			clprobdistPoissonInverseCDF(lambda, u, &err);
		}

		clock_t end = clock();
		float CPU_time = (float)(end - start) / CLOCKS_PER_SEC;
		printf("\n\tLambda : %1.1f (S) : Total CPU time (sec.) : %1.6f\n", lambda, CPU_time);
	}

}

void performanceOnGPU_Poisson(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, double lambda, size_t streamsBufSize, clrngMrg32k3aStream* streams, ExecType execType)
{
	printf("\tLambda : %1.1f (%s) : ", lambda, (execType == simWithObject ? "O" : "S"));
	cl_int err;
	size_t global_size = n1;

	//create stream buffer
	cl_mem stream_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams, &err);
	check_error(err, "%s(): cannot create stream buffer", __func__);

	clprobdistPoisson* poissonDist = NULL;
	cl_mem poissonDist_buffer = NULL;

	if (execType == simWithObject){
		//Create Poisson Dist
		clprobdistStatus err2;
		size_t poissonBufSize = 0;
		poissonDist = clprobdistPoissonCreate(lambda, &poissonBufSize, &err2);
		check_error(err, "%s(): cannot create Poisson Distribution", __func__);

		//Create Poisson buffer
		poissonDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, poissonBufSize, poissonDist, &err);
		check_error(err, "%s(): cannot create poissonDist_buffer buffer", __func__);

	}


	// Create the kernel.
	const char * preprocessors = "#define CLPROBDIST_NORMAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE";
	cl_program program = build_program_from_file(context, device,
		"client/Performance/PerformanceKernels.cl", preprocessors, PATH_RELATIVE_TO_LIB, BuildOption(_options, n / n1, lambda));


	cl_kernel kernel = clCreateKernel(program, kernelName, &err);
	check_error(err, "%s(): cannot create Kernel", __func__);

	// Set arguments for kernel and enqueue that kernel.
	size_t iarg = 0;
	err = clSetKernelArg(kernel, iarg++, sizeof(stream_buffer), &stream_buffer);
	if (execType == simWithObject){
		err |= clSetKernelArg(kernel, iarg++, sizeof(poissonDist_buffer), &poissonDist_buffer);
#ifdef USE_LOCAL_MEMORY
		err |= clSetKernelArg(kernel, iarg++, sizeof(poissonDist_buffer), NULL); //create local memory on the CU
#endif
	}
	//err |= clSetKernelArg(kernel, iarg++, sizeof(normalDist_buffer), &normalDist_buffer);
	check_error(err, "%s(): cannot create Kernel arguments", __func__);

	cl_event prof_event;        // Used for timing execution.

	// Enqueue kernels
	size_t max_wg_size = get_max_workgroup_size(device);
	size_t local_work_size = (global_size <= max_wg_size ? global_size : max_wg_size);

	err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, &local_work_size, 0, NULL, &prof_event);
	check_error(err, "%s(): cannot enqueue Kernel", __func__);

	// Finish processing the queue and compute the total time for executing the program.
	clFinish(queue);

	cl_ulong time_start, time_end;
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	cl_ulong total_GPU_time = time_end - time_start;

	printf("Total GPU time (sec.): %1.5f\n\n", (total_GPU_time / 1.0e9));


	// Deallocate resources
	clReleaseEvent(prof_event);
	clReleaseMemObject(stream_buffer);
	clReleaseMemObject(poissonDist_buffer);
	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clprobdistPoissonDestroy(poissonDist);

}


void performanceOnGPU_NoramlDist(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, size_t streamsBufSize, clrngMrg32k3aStream* streams, size_t normalBufSize, clprobdistNormal* normalDist)
{
	cl_int err;
	size_t global_size = n1;

	//create stream buffer
	cl_mem stream_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams, &err);
	check_error(err, "%s(): cannot create stream buffer", __func__);

	cl_mem normalDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, normalBufSize, normalDist, &err);
	check_error(err, "%s(): cannot create normalDist_buffer buffer", __func__);


	// Create the kernel.
	const char * preprocessors = "#define CLPROBDIST_NORMAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE";

	cl_program program = build_program_from_file(context, device,
		"client/Performance/PerformanceKernels.cl", preprocessors, PATH_RELATIVE_TO_LIB, BuildOption(_options, n / n1, 0));


	cl_kernel kernel = clCreateKernel(program, kernelName, &err);
	check_error(err, "%s(): cannot create Kernel", __func__);

	// Set arguments for kernel and enqueue that kernel.
	size_t iarg = 0;
	err = clSetKernelArg(kernel, iarg++, sizeof(stream_buffer), &stream_buffer);
	err |= clSetKernelArg(kernel, iarg++, sizeof(normalDist_buffer), &normalDist_buffer);

	check_error(err, "%s(): cannot create Kernel arguments", __func__);

	cl_event prof_event;        // Used for timing execution.

	// Enqueue kernels
	size_t max_wg_size = get_max_workgroup_size(device);
	size_t local_work_size = (global_size <= max_wg_size ? global_size : max_wg_size);

	err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, &local_work_size, 0, NULL, &prof_event);
	check_error(err, "%s(): cannot enqueue Kernel", __func__);

	// Finish processing the queue and compute the total time for executing the program.
	clFinish(queue);

	cl_ulong time_start, time_end;
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	cl_ulong total_GPU_time = time_end - time_start;

	printf("Total GPU time (sec.): %1.5f\n\n", (total_GPU_time / 1.0e9));


	// Deallocate resources
	clReleaseEvent(prof_event);
	clReleaseMemObject(stream_buffer);
	clReleaseMemObject(normalDist_buffer);
	clReleaseKernel(kernel);
	clReleaseProgram(program);

}

void performanceOnGPU_PoissonGlobalLocal(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	int n, int n1, double lambda, size_t streamsBufSize, clrngMrg32k3aStream* streams, ExecType execType)
{
	cl_int err;
	size_t global_size = n1;

	//create stream buffer
	cl_mem stream_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams, &err);
	check_error(err, "%s(): cannot create stream buffer", __func__);

	clprobdistPoisson* poissonDist = NULL;
	cl_mem poissonDist_buffer = NULL;


	//Create Poisson Dist
	clprobdistStatus err2;
	size_t poissonBufSize = 0;
	poissonDist = clprobdistPoissonCreate(lambda, &poissonBufSize, &err2);
	check_error(err, "%s(): cannot create Poisson Distribution", __func__);

	//Create Poisson buffer
	poissonDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, poissonBufSize, poissonDist, &err);
	check_error(err, "%s(): cannot create poissonDist_buffer buffer", __func__);

	// Create the kernel.
	cl_program program;
	char * preprocessors;
	const char * define1 = "#define CLPROBDIST_NORMAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE";
	const char * define2 = "#define CLPROBDIST_POISSON_OBJ_MEM CLPROBDIST_MEM_TYPE_LOCAL";

	if (execType == simUseLocalMem){
		preprocessors = (char *)malloc(strlen(define1) + strlen(define2) + 1);
		strcpy(preprocessors, define1);
		strcat(preprocessors, " \n");
		strcat(preprocessors, define2);

		program = build_program_from_file(context, device,
			"client/Performance/PerformanceKernels.cl", preprocessors, PATH_RELATIVE_TO_LIB, BuildOption(_options, n / n1, lambda));
		free(preprocessors);
	}
	else
		program = build_program_from_file(context, device,
			"client/Performance/PerformanceKernels.cl", define1, PATH_RELATIVE_TO_LIB, BuildOption(_options, n / n1, lambda));



	cl_kernel kernel = clCreateKernel(program, kernelName, &err);
	check_error(err, "%s(): cannot create Kernel", __func__);

	// Set arguments for kernel and enqueue that kernel.
	size_t iarg = 0;
	err = clSetKernelArg(kernel, iarg++, sizeof(stream_buffer), &stream_buffer);

	err |= clSetKernelArg(kernel, iarg++, sizeof(poissonDist_buffer), &poissonDist_buffer);
	if (execType == simUseLocalMem){
		err |= clSetKernelArg(kernel, iarg++, sizeof(poissonDist_buffer), NULL); //create local memory on the CU
	}
	check_error(err, "%s(): cannot create Kernel arguments", __func__);

	cl_event prof_event;        // Used for timing execution.

	// Enqueue kernels
	size_t max_wg_size = get_max_workgroup_size(device);
	size_t local_work_size = (global_size <= max_wg_size ? global_size : max_wg_size);

	err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_size, &local_work_size, 0, NULL, &prof_event);
	check_error(err, "%s(): cannot enqueue Kernel", __func__);

	// Finish processing the queue and compute the total time for executing the program.
	clFinish(queue);

	cl_ulong time_start, time_end;
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(prof_event, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	cl_ulong total_GPU_time = time_end - time_start;

	printf("Total GPU time (sec.): %1.5f\n\n", (total_GPU_time / 1.0e9));


	// Deallocate resources
	clReleaseEvent(prof_event);
	clReleaseMemObject(stream_buffer);
	clReleaseMemObject(poissonDist_buffer);
	clReleaseKernel(kernel);
	clReleaseProgram(program);
	clprobdistPoissonDestroy(poissonDist);

}
