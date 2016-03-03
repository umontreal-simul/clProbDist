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

#include "Types.h"
#include "../common.h"
#include <clRNG/mrg32k3a.h>

#include <clProbDist/clProbDist.h>
#include <clProbDist/poisson.h>
#include <clProbDist/normal.h>

/* [constants] */
#define L	100
#define c	2
#define h	0.1
#define K	10
#define k	1
#define per	0.95
/* [constants] *
 * If you uncomment the following, also uncomment
 * "#define USE_LOCAL_MEMORY" in InventoryKernels.cl.
 */
// #define USE_LOCAL_MEMORY

static char _options[255];

//************************************************************************
// Device Simulation
//************************************************************************

char* Build_Options(char _options[255], int m, int n2, int s, int S, int P, int selectedDist)
{
	char _m[7], _n2[7], _s[7], _S[7], _p[7], _selectedDist[7];
	sprintf(_m, "%d", m);
	sprintf(_n2, "%d", n2);
	sprintf(_s, "%d", s);
	sprintf(_S, "%d", S);
	sprintf(_p, "%d", P);
	sprintf(_selectedDist, "%d", selectedDist);

	cl_int err;
	strcpy(_options, clrngGetLibraryDeviceIncludes(&err));
	check_error(err, "cannot get clRNG device-side headers");

	strcat(_options, " -Dparam_m=\"");
	strcat(_options, _m);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_n2=\"");
	strcat(_options, _n2);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_s=\"");
	strcat(_options, _s);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_S=\"");
	strcat(_options, _S);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_P=\"");
	strcat(_options, _S);
	strcat(_options, "\"");

	strcat(_options, " -Dparam_SelectedDist=\"");
	strcat(_options, _selectedDist);

	return strcat(_options, "\"");
}

clrngStatus inventorySimulateRunsGPU(cl_context context, cl_device_id device, cl_command_queue queue,
	int m, int* s, int* S, int P, int n, int n1, int n2, ExecOption Option, const char * kernelName,
	size_t streamsBufSize, clrngMrg32k3aStream * streams_demand, clrngMrg32k3aStream * streams_order,
	int selectedDist, double* lambdaArr,
	size_t poissonBufSize, clprobdistPoisson* poissonDist, size_t normalBufSize, clprobdistNormal* normalDist,
	double *stat, simResult * results)
{
	cl_int err;
	size_t global_size;

	// Total number of work-items.
	if (n1 == 0) global_size = n;
	else{
		global_size = n1;
	}

	// Buffers to store the streams and profits.
	cl_mem streams_demand_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams_demand, &err);
	check_error(err, "%s(): cannot create streams demand buffer", __func__);

	cl_mem streams_order_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, streamsBufSize, streams_order, &err);
	check_error(err, "%s(): cannot create streams order buffer", __func__);

	cl_mem stat_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, n * sizeof(double), NULL, &err);
	check_error(err, "%s(): cannot create stat buffer", __func__);

	cl_mem poissonDist_buffer = NULL, normalDist_buffer = NULL, lambdaArr_buffer = NULL, s_buffer = NULL, S_buffer = NULL;
	if (selectedDist == 1){
		poissonDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, poissonBufSize, poissonDist, &err);
		check_error(err, "%s(): cannot create poissonDist_buffer buffer", __func__);

		lambdaArr_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
		normalDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
	}
	else if (selectedDist == 2){
		lambdaArr_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, m*sizeof(double), lambdaArr, &err);
		check_error(err, "%s(): cannot create lambdaArr_buffer buffer", __func__);

		poissonDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
		normalDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
	}
	else if (selectedDist == 3){
		normalDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, normalBufSize, normalDist, &err);
		check_error(err, "%s(): cannot create normalDist_buffer buffer", __func__);

		poissonDist_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
		lambdaArr_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(double), 0, &err);
	}

	//Buffers to store policies
	if (Option == Option2)
	{
		s_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, P * sizeof(int), s, &err);
		check_error(err, "%s(): cannot create streams demand buffer", __func__);

		S_buffer = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, P * sizeof(int), S, &err);
		check_error(err, "%s(): cannot create streams order buffer", __func__);
	}

	// Create kernel, that executes inventorySimulOne for each policy, on each work item.
	cl_program program = build_program_from_file(context, device,
		"client/Inventory/InventoryKernels.cl", NULL, PATH_RELATIVE_TO_LIB,
		Build_Options(_options, m, n2, s[0], S[0], P, selectedDist));

	cl_kernel kernel = clCreateKernel(program, kernelName, &err);
	check_error(err, "%s(): cannot create Kernel", __func__);

	// Set arguments for kernel and enqueue that kernel.
	size_t iarg = 0;
	err |= clSetKernelArg(kernel, iarg++, sizeof(streams_demand_buffer), &streams_demand_buffer);
	err |= clSetKernelArg(kernel, iarg++, sizeof(streams_order_buffer), &streams_order_buffer);
	err |= clSetKernelArg(kernel, iarg++, sizeof(stat_buffer), &stat_buffer);
	
	err |= clSetKernelArg(kernel, iarg++, sizeof(poissonDist_buffer), &poissonDist_buffer);
#ifdef USE_LOCAL_MEMORY
	err |= clSetKernelArg(kernel, iarg++, poissonBufSize, NULL); //create local memory on the CU
#endif
	err |= clSetKernelArg(kernel, iarg++, sizeof(lambdaArr_buffer), &lambdaArr_buffer);
#ifdef USE_LOCAL_MEMORY
	err |= clSetKernelArg(kernel, iarg++, m*sizeof(double), NULL); //create local memory on the CU
#endif
	err |= clSetKernelArg(kernel, iarg++, sizeof(normalDist_buffer), &normalDist_buffer);

	if (Option == Option2)
	{
		err |= clSetKernelArg(kernel, iarg++, sizeof(s_buffer), &s_buffer);
		err |= clSetKernelArg(kernel, iarg++, sizeof(S_buffer), &S_buffer);
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

	// After execution, transfer the results in the stat_diff array.
	err = clEnqueueReadBuffer(queue, stat_buffer, CL_TRUE, 0, n * sizeof(double), stat, 0, NULL, NULL);
	check_error(err, "%s(): cannot read the buffer", __func__);

	if ((Option != Option1) && (Option != Option2)) printf("\nTotal GPU time (sec.): %1.5f\n", (total_GPU_time / 1.0e9));
	printf("\n(s,S) = (%d,%d)\n", *s, *S);
	computeCI(n, stat, results);

	//Record results
	if (results != NULL){
		results->SimType = GPU_Exec;
		results->CPU_time = (total_GPU_time / 1.0e9);
	}

	//printf("***************************************************************************\n");

	// Deallocate resources
	clReleaseEvent(prof_event);
	clReleaseMemObject(stat_buffer);
	clReleaseMemObject(streams_demand_buffer);
	clReleaseMemObject(streams_order_buffer);
	clReleaseMemObject(s_buffer);
	clReleaseMemObject(S_buffer);
	clReleaseMemObject(poissonDist_buffer);
	clReleaseMemObject(lambdaArr_buffer);
	clReleaseMemObject(normalDist_buffer);
	clReleaseKernel(kernel);
	clReleaseProgram(program);

	return (clrngStatus)EXIT_SUCCESS;
}

//************************************************************************
// CPU Simulation
//************************************************************************
/* [simulate one run] */
double inventorySimulateOneRun(int m, int s, int S,
	double* lambdaArr, clrngMrg32k3aStream* stream_demand, clrngMrg32k3aStream*  stream_order,
	int selectedDist, clprobdistPoisson* poissonDist, clprobdistNormal* normalDist)
{
	clprobdistStatus err;
	double Xj = S, Yj;                // Stock in the morning and in the evening.
	double demand, profit = 0.0;     // Cumulated profit.
	double fa, fb, fbfa;

	if (selectedDist == 3){
		//Set Normal distribution trucated at [0, DBL_MAX[.
		fa = clprobdistNormalCDFWithObject(normalDist, 0, &err);
		fb = clprobdistNormalCDFWithObject(normalDist, DBL_MAX, &err);
		fbfa = fb - fa;
	}

	for (int j = 0; j < m; j++) {
		double u = clrngMrg32k3aRandomU01(stream_demand);
		//Poisson Dist using the object
		if (selectedDist == 1)
			demand =  clprobdistPoissonInverseCDFWithObject(poissonDist, u, &err);
			
		//Poisson Dist using the static method
		else if (selectedDist == 2)
			demand = clprobdistPoissonInverseCDF(lambdaArr[j], u, &err);
	
		//Normal Dist using the object
		else if (selectedDist == 3)
			demand = clprobdistNormalInverseCDFWithObject(normalDist, fa + fbfa * u, &err) * clprobdistNormalGetSigma(normalDist, &err) + clprobdistNormalGetMu(normalDist, &err);

		Yj = Xj - demand;
		if (Yj < 0)
			Yj = 0;                  // Lost demand.
		profit += c * (Xj - Yj) - h * Yj;
		double prob = clrngMrg32k3aRandomU01(stream_order);

		if ((Yj < s) && (prob < per)) {
			// We have a successful order.
			profit -= K + k * (S - Yj);
			Xj = S;
		}
		else
			Xj = Yj;

	}

	return profit / m;
}
/* [simulate one run] */

clrngStatus inventorySimulateRunsCPU(int m, int s, int S, int n,
	double* lambdaArr, clrngMrg32k3aStream* stream_demand, clrngMrg32k3aStream* stream_order,
	int selectedDist, clprobdistPoisson* poissonDist, clprobdistNormal* normalDist,
	double *stat_profit, ExecType execType, simResult * results)
{
	clock_t start = clock();

	if (execType == basic){
		// (basic case) : Performs n independent simulation runs of the system for m days with the (s,S) policy
		// using a single stream and the same substream for everything, and saves daily profit values.
		// Equivalent implementation of inventorySimulateRunsOneStream() from the document.
		for (int i = 0; i < n; i++)
			stat_profit[i] = inventorySimulateOneRun(m, s, S,
			lambdaArr, stream_demand, stream_demand,
			selectedDist, poissonDist, normalDist);
		if (results != NULL) results->ExecOption = 1;
	}

	else if (execType == Case_a){
		//Case (a) : Similar to (basic), but using two streams and their substreams.
		// Equivalent implementation of inventorySimulateRunsSubstreams() from the document.
		for (int i = 0; i < n; i++){
			stat_profit[i] = inventorySimulateOneRun(m, s, S, lambdaArr,
				stream_demand, stream_order,
				selectedDist, poissonDist, normalDist);
			clrngMrg32k3aForwardToNextSubstreams(1, stream_demand);
			clrngMrg32k3aForwardToNextSubstreams(1, stream_order);
		}
		if (results != NULL) results->ExecOption = 2;
	}

	else if (execType == Case_b){
		//Case (b) : Similar to (basic), but with two arrays of n streams each, using a new pair of streams for each run.
		// Equivalent implementation of inventorySimulateRunsManyStreams() from the document.
		for (int i = 0; i < n; i++)
			stat_profit[i] = inventorySimulateOneRun(m, s, S,
			lambdaArr, &stream_demand[i], &stream_order[i],
			selectedDist, poissonDist, normalDist);
		if (results != NULL) results->ExecOption = 3;
	}
	//Compute Execution Time
	clock_t end = clock();
	float CPU_time = (float)(end - start) / CLOCKS_PER_SEC;
	printf("\nTotal CPU time (sec.): %1.6f\n", CPU_time);

	//Compute CI
	computeCI(n, stat_profit, results);

	//Record results
	if (results != NULL){
		results->SimType = CPU_Exec;
		results->CPU_time = CPU_time;
	}

	return (clrngStatus)EXIT_SUCCESS;
}
