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
#include <math.h>
#include <string.h>
#include <clRNG/clRNG.h>
#include "Types.h"
#include "../common.h"

#include <clRNG/mrg32k3a.h>
#include <clProbDist/normal.h>

static char _options[255];
char* BuildOptions(char _options[255], int n2, int s, double strike, double discount)
{
	char _n2[7], _s[7], _strike[20], _discount[20];

	sprintf(_n2, "%d", n2);
	sprintf(_s, "%d", s);
	sprintf(_strike, "%f", strike);
	sprintf(_discount, "%f", discount);

	cl_int err;
	strcpy(_options, clrngGetLibraryDeviceIncludes(&err));
	check_error(err, "cannot get clRNG device-side headers");

	strcat(_options, " -Dparam_n2=\"");
	strcat(_options, _n2);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_s=\"");
	strcat(_options, _s);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_strike=\"");
	strcat(_options, _strike);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_discount=\"");
	strcat(_options, _discount);

	return strcat(_options, "\"");
}

int AsianOptionSimulateRunsGPU(cl_context context, cl_device_id device, cl_command_queue queue, const char * kernelName,
	                           int n, int n1, int s, double strike, double discount, double* logS, double* muDelta, double* sigmaSqrtDelta,
							   size_t streamsBufSize, clrngMrg32k3aStream * streams,
							   size_t normalBufSize, clprobdistNormal* normalDist, double *stat)

{
	cl_int err;

	// Total number of work-items.
	size_t global_size = n1;
	int n2 = n / n1;

	// Buffers to store the streams and profits.
	cl_mem streams_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams, &err);
	check_error(err, "%s(): cannot create streams buffer", __func__);

	cl_mem normalDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, normalBufSize, normalDist, &err);
	check_error(err, "%s(): cannot create normalDist_buffer buffer", __func__);

	cl_mem stat_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, n * sizeof(double), NULL, &err);
	check_error(err, "%s(): cannot create stat buffer", __func__);

	cl_mem logS_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (s + 1) * sizeof(double), logS, &err);
	check_error(err, "%s(): cannot create logS buffer", __func__);

	cl_mem muDelta_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, s * sizeof(double), muDelta, &err);
	check_error(err, "%s(): cannot create muDelta buffer", __func__);

	cl_mem sigmaSqrtDelta_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, s * sizeof(double), sigmaSqrtDelta, &err);
	check_error(err, "%s(): cannot create sigmaSqrtDelta buffer", __func__);

	// Create kernel, that executes inventorySimulOne for each policy, on each work item.
	cl_program program = build_program_from_file(context, device,
		"client/AsianOption/AsianOptionKernel.cl", NULL, PATH_RELATIVE_TO_LIB,
		BuildOptions(_options, n2, s, strike, discount));

	cl_kernel kernel = clCreateKernel(program, kernelName, &err);
	check_error(err, "%s(): cannot create Kernel", __func__);

	// Set arguments for kernel and enqueue that kernel.
	err |= clSetKernelArg(kernel, 0, sizeof(streams_buffer), &streams_buffer);
	err |= clSetKernelArg(kernel, 1, sizeof(normalDist_buffer), &normalDist_buffer);
	err |= clSetKernelArg(kernel, 2, sizeof(stat_buffer), &stat_buffer);

	err |= clSetKernelArg(kernel, 3, sizeof(logS_buffer), &logS_buffer);
	err |= clSetKernelArg(kernel, 4, sizeof(muDelta_buffer), &muDelta_buffer);
	err |= clSetKernelArg(kernel, 5, sizeof(sigmaSqrtDelta_buffer), &sigmaSqrtDelta_buffer);

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

	// After execution, transfer the results in the stat array.
	err = clEnqueueReadBuffer(queue, stat_buffer, CL_TRUE, 0, n * sizeof(double), stat, 0, NULL, NULL);
	check_error(err, "%s(): cannot read the buffer", __func__);

	printf("\nTotal GPU time (sec.): %1.5f\n", (total_GPU_time / 1.0e9));
	computeCI(n, stat, NULL);

	// Deallocate resources
	clReleaseEvent(prof_event);
	clReleaseMemObject(streams_buffer);
	clReleaseMemObject(normalDist_buffer);
	clReleaseMemObject(stat_buffer);
	clReleaseMemObject(logS_buffer);
	clReleaseMemObject(muDelta_buffer);
	clReleaseMemObject(sigmaSqrtDelta_buffer);
	clReleaseKernel(kernel);
	clReleaseProgram(program);

	return EXIT_SUCCESS;
}

double GetFloatPrecision(double x, double precision)
{
	return (floor((x * pow(10, precision) + 0.5)) / pow(10, precision));
}

// Generates the process S.
void generatePath(int s, double* logS, double* muDelta, double* sigmaSqrtDelta, clrngMrg32k3aStream* stream, clprobdistNormal* normalDist) {
	clprobdistStatus err;

	for (int j = 0; j < s; j++){

		double u = clrngMrg32k3aRandomU01(stream);
		double d = clprobdistNormalInverseCDFWithObject(normalDist, u, &err);

		logS[j + 1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j] * d;
	}

}

// Computes and returns the discounted option payoff.
double getPayoff(int s, double strike, double discount, double* logS) {
	double average = 0.0;  // Average of the GBM process.
	for (int j = 1; j <= s; j++) average += exp(logS[j]);
	average /= s;
	if (average > strike) return discount * (average - strike);
	else return 0.0;
}

// Performs n indep. run
int SimulateAsianOption(cl_context context, cl_device_id device, cl_command_queue queue, void* data_)
{
	clrngStatus err;
	clprobdistStatus err2;

	AsianOptionData* data = (AsianOptionData*)data_;
	int n = data->n;
	int n1 = data->n1;
	double r = data->r;
	double sigma = data->sigma;
	double strike = data->strike;
	double s0 = data->s0;
	int s = data->s;
	double* zeta = data->zeta;

	//Set Simulation params
	double discount = GetFloatPrecision(exp(-r * zeta[s]), 2);

	double mu = r - 0.5 * sigma * sigma;

	double* muDelta = (double*)malloc(s*sizeof(double));
	double* sigmaSqrtDelta = (double*)malloc(s*sizeof(double));
	double* logS = (double*)malloc((s + 1)*sizeof(double));

	double delta;
	for (int j = 0; j < s; j++) {
		delta = zeta[j + 1] - zeta[j];
		muDelta[j] = mu * delta;
		sigmaSqrtDelta[j] = sigma * sqrt(delta);
	}
	logS[0] = log(s0);

	//Create payOff stat
	double *stat_payOff = (double *)malloc(n * sizeof(double));

	//Create Stream 
	size_t streamBufferSize = 0;
	clrngMrg32k3aStream* streams = clrngMrg32k3aCreateStreams(NULL, n1, &streamBufferSize, &err);
	check_error(err, "%s(): cannot create random stream", __func__);

	//Set Std Noraml Distribution
	size_t normalBufSize = 0;
	clprobdistNormal* normalDist = clprobdistNormalCreate(0, 1, &normalBufSize, &err2);


	//*******************************************************************************
	//Simulate on CPU
	//*******************************************************************************
	int n2 = n / n1;
	clock_t start = clock();
	for (int j = 0; j < n1; j++) {
		clrngMrg32k3aStream* stream = &streams[j];

		for (int i = 0; i < n2; i++) {
			generatePath(s, logS, muDelta, sigmaSqrtDelta, stream, normalDist);
			stat_payOff[i*n1 + j] = getPayoff(s, strike, discount, logS);
			
			clrngMrg32k3aForwardToNextSubstreams(1, stream);
		}
	}

	//Compute Execution Time
	clock_t end = clock();
	float CPU_time = (float)(end - start) / CLOCKS_PER_SEC;
	printf("\nTotal CPU time (sec.): %1.6f\n", CPU_time);

	//Compute CI
	computeCI(n, stat_payOff, NULL);

	//*******************************************************************************
	//Simulate on GPU
	//*******************************************************************************
	clrngMrg32k3aRewindStreams(n1, streams);

	AsianOptionSimulateRunsGPU(context, device, queue, "AsianOptionSimulateGPU",
                         	   n, n1, s, strike, discount, logS, muDelta, sigmaSqrtDelta,
							   streamBufferSize, streams,
		                       normalBufSize, normalDist,stat_payOff);

	//Free resources
	free(stat_payOff);
	free(muDelta);
	free(sigmaSqrtDelta);
	free(logS);

	clrngMrg32k3aDestroyStreams(streams);
	clprobdistNormalDestroy(normalDist);

	return EXIT_SUCCESS;

}

