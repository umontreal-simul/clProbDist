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

#define CLPROBDIST_NORMAL_OBJ_MEM CLPROBDIST_MEM_TYPE_PRIVATE
#define CLRNG_ENABLE_SUBSTREAMS

/*! [clRNG header] */
#include <clRNG/mrg32k3a.clh>
/*! [clRNG header] */

#include "clProbDist/normal.clh"

#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

// Generates the process S.
void generatePath(int s, __global double* logS, __global double* muDelta, __global double* sigmaSqrtDelta, 
	              clrngMrg32k3aStream* stream, clprobdistNormal* normalDist) {
	clprobdistStatus err;

	for (int j = 0; j < s; j++){

		double u = clrngMrg32k3aRandomU01(stream);
		double d = clprobdistNormalInverseCDFWithObject(normalDist, u, &err);

		logS[j + 1] = logS[j] + muDelta[j] + sigmaSqrtDelta[j] * d;
	}

}

// Computes and returns the discounted option payoff.
double getPayoff(int s, double strike, double discount, __global double* logS) {
	double average = 0.0;  // Average of the GBM process.
	for (int j = 1; j <= s; j++) average += exp(logS[j]);
	average /= s;
	if (average > strike) return discount * (average - strike);
	else return 0.0;
}

//use distinct streams across the work items and n2 substreams within each work item.
__kernel void AsianOptionSimulateGPU(__global clrngMrg32k3aHostStream* streams, __global clprobdistNormal* g_normalDist, __global double* stat_payOff,
	                                 __global double* logS, __global double* muDelta, __global double* sigmaSqrtDelta)
{
	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.
	int n1 = get_global_size(0); //Total number of work items

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

	//Make a copy of clprobdistNormal in private memory
	clprobdistNormal p_normalDist = *g_normalDist;

	//Simulate
	for (int i = 0; i < param_n2; i++) {

		generatePath(param_s, logS, muDelta, sigmaSqrtDelta, &stream_d, &p_normalDist);
		stat_payOff[i * n1 + gid] = getPayoff(param_s, param_strike, param_discount, logS);

		clrngMrg32k3aForwardToNextSubstreams(1, &stream_d);	
	}
}
