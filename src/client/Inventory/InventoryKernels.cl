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

#define CLRNG_ENABLE_SUBSTREAMS
#define CLPROBDIST_NORMAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE

//Important : all the variables that start with the prefix 'param_' are passed to the Open CL C compiler as inline parameter (see.: BuildOptions)

/*
 * If you uncomment the following, also uncomment
 * "#define USE_LOCAL_MEMORY" in SimulateRuns.c.
 */
// #define CLPROBDIST_POISSON_OBJ_MEM CLPROBDIST_MEM_TYPE_LOCAL

/*! [clRNG header] */
#include <clRNG/mrg32k3a.clh>
/*! [clRNG header] */

#include "clProbDist/poisson.clh"
#include "clProbDist/normal.clh"

__constant int L = 100;
__constant int c = 2;
__constant double h = 0.1;
__constant int K = 10;
__constant int k = 1;
__constant double p = 0.95;


#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
clprobdistStatus CopyLambdaArrayFromGlobal(cl_uint len, __local double* dest, __global double* src)
{
	//Check params
	if (!dest)
		return clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "CopyLambdaArrayFromGlobal(): dest cannot be NULL");
	if (!src)
		return clprobdistSetErrorString(CLPROBDIST_INVALID_VALUE, "CopyLambdaArrayFromGlobal(): src cannot be NULL");

	//Copy values from global to local memory
	for (size_t i = 0; i < len; i++)
		dest[i] = src[i];

	return CLPROBDIST_SUCCESS;
}
#endif



double inventorySimulateOneRunFixedPoissonLambda(int m, int s, int S,
	clrngMrg32k3aStream* stream_demand, clrngMrg32k3aStream* stream_order,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__local clprobdistPoisson* poissonDist
#else
	__global clprobdistPoisson* poissonDist
#endif
	)
{
	clprobdistStatus err;
	double Xj = S, Yj;                // Stock in the morning and in the evening.
	double demand, profit = 0.0;     // Cumulated profit.

	for (int j = 0; j < m; j++) {
		// Generate demand
		cl_double u = clrngMrg32k3aRandomU01(stream_demand);

		//Call to method that uses the Poisson distribution object
		demand = clprobdistPoissonInverseCDFWithObject(poissonDist, u, &err);

		// Subtract the demand for the day.		 
		Yj = Xj - demand;
		if (Yj < 0)
			Yj = 0;                  // Lost demand.
		profit += c * (Xj - Yj) - h * Yj;

		double prob = clrngMrg32k3aRandomU01(stream_order);
		//printf("%f \n", prob);

		if ((Yj < s) && (prob < p)) {
			// We have a successful order.
			profit -= K + k * (S - Yj);
			Xj = S;
		}
		else
			Xj = Yj;
	}
	return profit / m;

}

double inventorySimulateOneRunDynamicLambdaPoisson(int m, int s, int S,
	clrngMrg32k3aStream* stream_demand, clrngMrg32k3aStream* stream_order,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__local double* lambdaArray
#else
	__global double* lambdaArray
#endif
	)
{
	clprobdistStatus err;
	double Xj = S, Yj;                // Stock in the morning and in the evening.
	double demand, profit = 0.0;     // Cumulated profit.

	for (int j = 0; j < m; j++) {
		// Generate demand
		cl_double u = clrngMrg32k3aRandomU01(stream_demand);

		//Call to static version of InverseCDF function
		demand = clprobdistPoissonInverseCDF(lambdaArray[j], u, &err);

		// Subtract the demand for the day.		 
		Yj = Xj - demand;
		if (Yj < 0)
			Yj = 0;                  // Lost demand.
		profit += c * (Xj - Yj) - h * Yj;

		double prob = clrngMrg32k3aRandomU01(stream_order);

		if ((Yj < s) && (prob < p)) {
			// We have a successful order.
			profit -= K + k * (S - Yj);
			Xj = S;
		}
		else
			Xj = Yj;
	}
	return profit / m;

}

double inventorySimulateOneRunNormal(int m, int s, int S,
	clrngMrg32k3aStream* stream_demand, clrngMrg32k3aStream* stream_order,
	clprobdistNormal* normalDist)
{
	clprobdistStatus err;
	double Xj = S, Yj;                // Stock in the morning and in the evening.
	double demand, profit = 0.0;     // Cumulated profit.

	//Set Normal distribution trucated at [0, DBL_MAX[.
	double fa = clprobdistNormalCDFWithObject(normalDist, 0, &err);
	double fb = clprobdistNormalCDFWithObject(normalDist, DBL_MAX, &err);
	double fbfa = fb - fa;


	for (int j = 0; j < m; j++) {
		// Generate demand
		cl_double u = clrngMrg32k3aRandomU01(stream_demand);

		demand = clprobdistNormalInverseCDFWithObject(normalDist, fa + fbfa * u, &err) * clprobdistNormalGetSigma(normalDist, &err) + clprobdistNormalGetMu(normalDist, &err);

		// Subtract the demand for the day.		 
		Yj = Xj - demand;
		if (Yj < 0)
			Yj = 0;                  // Lost demand.
		profit += c * (Xj - Yj) - h * Yj;

		double prob = clrngMrg32k3aRandomU01(stream_order);
		//printf("%f \n", prob);

		if ((Yj < s) && (prob < p)) {
			// We have a successful order.
			profit -= K + k * (S - Yj);
			Xj = S;
		}
		else
			Xj = Yj;
	}
	return profit / m;

}

//************************************************************************
// One policy
//************************************************************************

//Case (a) and (b) : Simulate n runs on n work items using two streams and their substreams or 2*n streams
__kernel void inventorySimulateGPU(__global clrngMrg32k3aHostStream* streams_demand, __global clrngMrg32k3aHostStream* streams_order, __global double* stat_profit,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__global clprobdistPoisson* g_poissonDist,
	__local  clprobdistPoisson* poissonDist,
	__global double* g_lambdaArr,
	__local  double* lambdaArr,
#else
	__global clprobdistPoisson* poissonDist,
	__global double* lambdaArr,
#endif
	__global clprobdistNormal* g_normalDist
	)
{
	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.
	int group_size = get_local_size(0);
	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_demand_d, stream_order_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[gid]);
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[gid]);

	clprobdistNormal p_normalDist;
	if (param_SelectedDist == 1){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global Poissont Distribution to local memory.
		if (gid % group_size == 0)
			clprobdistPoissonCopyOverFromGlobal(poissonDist, g_poissonDist);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		stat_profit[gid] = inventorySimulateOneRunFixedPoissonLambda(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, poissonDist);
		 
	}
	else if (param_SelectedDist == 2){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global LambdArray to local memory.
		if (gid % group_size == 0)
			CopyLambdaArrayFromGlobal(param_m, lambdaArr, g_lambdaArr);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		stat_profit[gid] = inventorySimulateOneRunDynamicLambdaPoisson(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, lambdaArr);
	}
	else if (param_SelectedDist == 3){
          p_normalDist = *g_normalDist;

		  stat_profit[gid] = inventorySimulateOneRunNormal(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, &p_normalDist);
	}
}

//Case (c) : use distinct streams across the work items and n2 substreams within each work item.
__kernel void inventorySimulSubstreamsGPU(__global clrngMrg32k3aHostStream *streams_demand, __global clrngMrg32k3aHostStream *streams_order, __global double *stat_profit,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__global clprobdistPoisson* g_poissonDist,
	__local  clprobdistPoisson* poissonDist,
	__global double* g_lambdaArr,
	__local  double* lambdaArr,
#else
	__global clprobdistPoisson* poissonDist,
	__global double* lambdaArr,
#endif
	__global clprobdistNormal* g_normalDist
	)
{
	// Each of the n1 work items executes the following to simulate n2 runs.
	int gid = get_global_id(0); // Id of this work item
	int n1 = get_global_size(0); //Total number of work items
	int group_size = get_local_size(0);

	// Make local copies of the stream states in private memory
	clrngMrg32k3aStream stream_demand_d, stream_order_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[gid]);
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[gid]);

	clprobdistNormal p_normalDist;
	if (param_SelectedDist == 1){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global Poissont Distribution to local memory.
		if (gid % group_size == 0)
			clprobdistPoissonCopyOverFromGlobal(poissonDist, g_poissonDist);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++) {
			stat_profit[i * n1 + gid] = inventorySimulateOneRunFixedPoissonLambda(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, poissonDist);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
	else if (param_SelectedDist == 2){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global LambdArray to local memory.
		if (gid % group_size == 0)
			CopyLambdaArrayFromGlobal(param_m, lambdaArr, g_lambdaArr);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++) {
			stat_profit[i * n1 + gid] = inventorySimulateOneRunDynamicLambdaPoisson(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, lambdaArr);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
	else if (param_SelectedDist == 3){
        p_normalDist = *g_normalDist;

		for (int i = 0; i < param_n2; i++) {
			stat_profit[i * n1 + gid] = inventorySimulateOneRunNormal(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, &p_normalDist);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
}

//Case (d) : each work item uses 2*n2 distinct streams instead of using substreams.
__kernel void inventorySimul_DistinctStreams_GPU(__global clrngMrg32k3aHostStream *streams_demand, __global clrngMrg32k3aHostStream *streams_order, __global double *stat_profit,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__global clprobdistPoisson* g_poissonDist,
	__local  clprobdistPoisson* poissonDist,
	__global double* g_lambdaArr,
	__local  double* lambdaArr,
#else
	__global clprobdistPoisson* poissonDist,
	__global double* lambdaArr,
#endif
	__global clprobdistNormal* g_normalDist
	)
{
	// Each of the n1 work items executes the following to simulate n2 runs.
	int gid = get_global_id(0); // Id of this work item.
	int n1 = get_global_size(0); //Total number of work items
	int group_size = get_local_size(0);

	clrngMrg32k3aStream stream_demand_d, stream_order_d;

	clprobdistNormal p_normalDist;
	if (param_SelectedDist == 1){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global Poissont Distribution to local memory.
		if (gid % group_size == 0)
			clprobdistPoissonCopyOverFromGlobal(poissonDist, g_poissonDist);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++) {
			// Make local copies of the stream states, in private memory.
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[i * n1 + gid]);
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[i * n1 + gid]);

			stat_profit[i * n1 + gid] = inventorySimulateOneRunFixedPoissonLambda(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, poissonDist);
		}
	}
	else if (param_SelectedDist == 2){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global LambdArray to local memory.
		if (gid % group_size == 0)
			CopyLambdaArrayFromGlobal(param_m, lambdaArr, g_lambdaArr);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++) {
			// Make local copies of the stream states, in private memory.
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[i * n1 + gid]);
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[i * n1 + gid]);

			stat_profit[i * n1 + gid] = inventorySimulateOneRunDynamicLambdaPoisson(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, lambdaArr);
		}
	}
	else if (param_SelectedDist == 3){
        p_normalDist = *g_normalDist;
		for (int i = 0; i < param_n2; i++) {
			// Make local copies of the stream states, in private memory.
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[i * n1 + gid]);
			clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[i * n1 + gid]);

			stat_profit[i * n1 + gid] = inventorySimulateOneRunNormal(param_m, param_s, param_S, &stream_demand_d, &stream_order_d, &p_normalDist);
		}
	}

}

//************************************************************************
// Several policies
//************************************************************************

//Simulate n2 runs on n1p workitmes using n1 streams and their substreams
__kernel void inventorySimulPoliciesGPU(__global clrngMrg32k3aHostStream *streams_demand, __global clrngMrg32k3aHostStream *streams_order, __global double *stat_profit,
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__global clprobdistPoisson* g_poissonDist,
	__local  clprobdistPoisson* poissonDist,
	__global double* g_lambdaArr,
	__local  double* lambdaArr,
#else
	__global clprobdistPoisson* poissonDist,
	__global double* lambdaArr,
#endif
	__global clprobdistNormal* g_normalDist,
	__global cl_int *s, __global cl_int *S)
{
	// Each of the n1*P work items executes the following to simulate n2 runs.
	int gid = get_global_id(0);    // Id of this work item.
	int n1p = get_global_size(0);  // n1*P = Total number of work items
	int n1 = n1p / param_P;       // Number of streams.
	int k = gid / n1;            // Index that identify which policy the work item will use. in case p = 2 : k = 0 or k = 1
	int j = gid % n1;            // Index that identify the index of the work item modulo n1, in case n1=100 : j=0..99
	int group_size = get_local_size(0);

	// Make local copies of the stream states, in private memory.
	clrngMrg32k3aStream stream_demand_d, stream_order_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_demand_d, &streams_demand[j]);
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_order_d, &streams_order[j]);

	clprobdistNormal p_normalDist;
	if (param_SelectedDist == 1){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global Poissont Distribution to local memory.
		if (gid % group_size == 0)
			clprobdistPoissonCopyOverFromGlobal(poissonDist, g_poissonDist);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++)
		{
			stat_profit[i * n1p + gid] = inventorySimulateOneRunFixedPoissonLambda(param_m, s[k], S[k], &stream_demand_d, &stream_order_d, poissonDist);

			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
	else if (param_SelectedDist == 2){
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
		//First work-item of each group make a copy of global LambdArray to local memory.
		if (gid % group_size == 0)
			CopyLambdaArrayFromGlobal(param_m, lambdaArr, g_lambdaArr);

		//All work-items sync here to make sure that the first work-item has finished copying.
		barrier(CLK_LOCAL_MEM_FENCE);
#endif

		for (int i = 0; i < param_n2; i++)
		{
			stat_profit[i * n1p + gid] = inventorySimulateOneRunDynamicLambdaPoisson(param_m, s[k], S[k], &stream_demand_d, &stream_order_d, lambdaArr);

			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
	else if (param_SelectedDist == 3){
		p_normalDist = *g_normalDist;

		for (int i = 0; i < param_n2; i++)
		{
			stat_profit[i * n1p + gid] = inventorySimulateOneRunNormal(param_m, s[k], S[k], &stream_demand_d, &stream_order_d, &p_normalDist);

			clrngMrg32k3aForwardToNextSubstreams(1, &stream_demand_d);
			clrngMrg32k3aForwardToNextSubstreams(1, &stream_order_d);
		}
	}
}

