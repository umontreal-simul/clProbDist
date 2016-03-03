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

/*! [clRNG header] */
#include <clRNG/mrg32k3a.clh>
/*! [clRNG header] */

#include "clProbDist/poisson.clh"
#include "clProbDist/normal.clh"

#pragma OPENCL EXTENSION cl_amd_printf : enable


//******************************************************************************************
// Option 1 :
//******************************************************************************************
//use distinct streams across the work items and n2 loop within each work item.
__kernel void performancePoissonWithObject(__global clrngMrg32k3aHostStream* streams,_CLPROBDIST_POISSON_OBJ_MEM clprobdistPoisson* poissonDist)
{
	clprobdistStatus err;

	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.
	int group_size = get_local_size(0);

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

	//Simulate n2 random nbr generation
	for (int i = 0; i < param_n2; i++) {
		clprobdistPoissonInverseCDFWithObject(poissonDist, clrngMrg32k3aRandomU01(&stream_d), &err);
	}
}

//use distinct streams across the work items and n2 loop within each work item.
__kernel void performancePoisson(__global clrngMrg32k3aHostStream* streams)
{
	clprobdistStatus err;

	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

	//Simulate n2 random nbr generation
	for (int i = 0; i < param_n2; i++) {
		clprobdistPoissonInverseCDF(param_lambda, clrngMrg32k3aRandomU01(&stream_d), &err);
	}
}

//******************************************************************************************
// Option 2 :
//******************************************************************************************

__kernel void performancePrivateNormalDist(__global clrngMrg32k3aHostStream* streams, __global clprobdistNormal* g_normalDist)
{
	clprobdistStatus err;

	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

	//Make a copy of clprobdistNormal in private memory
	clprobdistNormal p_normalDist = *g_normalDist;

	//Simulate n2 random nbr generation
	for (int i = 0; i < param_n2; i++) {
		clprobdistNormalInverseCDFWithObject(&p_normalDist, clrngMrg32k3aRandomU01(&stream_d), &err);
	}
}
__kernel void performanceGlobalNormalDist(__global clrngMrg32k3aHostStream* streams, __constant clprobdistNormal* c_normalDist)
{
	clprobdistStatus err;

	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

	//Simulate n2 random nbr generation
	for (int i = 0; i < param_n2; i++) {
		clprobdistNormalInverseCDFWithObject(&c_normalDist, clrngMrg32k3aRandomU01(&stream_d), &err);
	}
}

//******************************************************************************************
// Option 3 :
//******************************************************************************************

__kernel void performancePoissonGlobalLocal(__global clrngMrg32k3aHostStream* streams, 
#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	__global clprobdistPoisson* g_poissonDist,
	__local  clprobdistPoisson* poissonDist
#else
	_CLPROBDIST_POISSON_OBJ_MEM clprobdistPoisson* poissonDist
#endif
)
{
	clprobdistStatus err;

	// Each of the n work items executes the following code.
	int gid = get_global_id(0); // Id of this work item.
	int group_size = get_local_size(0);

	// Make copies of the stream states in private memory.
	clrngMrg32k3aStream stream_d;
	clrngMrg32k3aCopyOverStreamsFromGlobal(1, &stream_d, &streams[gid]);

#if CLPROBDIST_POISSON_OBJ_MEM == CLPROBDIST_MEM_TYPE_LOCAL
	//First work-item of each group make a copy of global Poissont Distribution to local memory.
	if (gid % group_size == 0)
		clprobdistPoissonCopyOverFromGlobal(poissonDist, g_poissonDist);

	//All work-items sync here to make sure that the first work-item has finished copying.
	barrier(CLK_LOCAL_MEM_FENCE);
#endif

	//Simulate n2 random nbr generation
	for (int i = 0; i < param_n2; i++) {
		clprobdistPoissonInverseCDFWithObject(poissonDist, clrngMrg32k3aRandomU01(&stream_d), &err);
	}
}
