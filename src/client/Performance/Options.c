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
#include <clRNG/clRNG.h>
#include "../common.h"
#include "Options.h"
#include "SimulateRun.h"

/*! [clRNG header] */
#include <clRNG/mrg31k3p.h>
/*! [clRNG header] */

#include "../../include/clProbDist/gamma.h"
#include "../../include/clProbDist/normal.h"

int Option1(cl_context context, cl_device_id device, cl_command_queue queue, void* data_)
{
	OptionData* data = (OptionData*)data_;
	int n = data->n;
	int n1 = data->n1;
	double * lambdaArr = data->lambdaArr;

	//Declare vars 
	clrngStatus err;
	size_t streamBufferSize;

	//Create stream demand
	clrngMrg32k3aStream* streams = clrngMrg32k3aCreateStreams(NULL, n1, &streamBufferSize, &err);
	check_error(err, "%s(): cannot create random streams demand", __func__);

	printf("\nSimulate Poisson Distribution with Object :\n");
	printf("On Device :\n");
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			performanceOnGPU_Poisson(context, device, queue, "performancePoissonWithObject", n, n1, lambdaArr[i], streamBufferSize, streams, simWithObject);

	//**************************************************************
	printf("\nSimulate Poisson Distribution with static calls :\n");
	printf("On Device :\n");
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			performanceOnGPU_Poisson(context, device, queue, "performancePoisson", n, n1, lambdaArr[i], streamBufferSize, streams, simStatic);
	
	//Free resources
	clrngMrg32k3aDestroyStreams(streams);

	return EXIT_SUCCESS;
}

int Option2(cl_context context, cl_device_id device, cl_command_queue queue, void* data_)
{
	OptionData* data = (OptionData*)data_;
	int n = data->n;
	int n1 = data->n1;

	//Declare vars 
	clrngStatus err;
	size_t streamBufferSize;

	//Create stream demand
	clrngMrg32k3aStream* streams = clrngMrg32k3aCreateStreams(NULL, n1, &streamBufferSize, &err);
	check_error(err, "%s(): cannot create random streams demand", __func__);

	//Set Std Noraml Distribution
	size_t normalBufSize = 0;
	clprobdistNormal* normalDist = clprobdistNormalCreate(0, 1, &normalBufSize, NULL);

	printf("\nSimulate Normal Distribution copied on private memory :\n");
	for (int i = 0; i < 3; i++)
		performanceOnGPU_NoramlDist(context, device, queue, "performancePrivateNormalDist", n, n1, streamBufferSize, streams, normalBufSize, normalDist);

	//**************************************************************
	printf("\n\nSimulate Normal Distribution copied on global memory :\n");
	printf("On Device :\n");
	for (int i = 0; i < 3; i++)
		performanceOnGPU_NoramlDist(context, device, queue, "performanceGlobalNormalDist", n, n1, streamBufferSize, streams, normalBufSize, normalDist);

	//Free resources
	clrngMrg32k3aDestroyStreams(streams);
	clprobdistNormalDestroy(normalDist);
	return EXIT_SUCCESS;
}

int Option3(cl_context context, cl_device_id device, cl_command_queue queue, void* data_)
{
	OptionData* data = (OptionData*)data_;
	int n = data->n;
	int n1 = data->n1;
	double * lambdaArr = data->lambdaArr;

	//Declare vars 
	clrngStatus err;
	size_t streamBufferSize;

	//Create stream demand
	clrngMrg32k3aStream* streams = clrngMrg32k3aCreateStreams(NULL, n1, &streamBufferSize, &err);
	check_error(err, "%s(): cannot create random streams demand", __func__);

	printf("\nSimulate Poisson Distribution with Object in local memory :\n");
	for (int i = 0; i < 3; i++)
		performanceOnGPU_PoissonGlobalLocal(context, device, queue, "performancePoissonGlobalLocal", n, n1, lambdaArr[0], streamBufferSize, streams, simUseLocalMem);

	//**************************************************************
	printf("\nSimulate Poisson Distribution with Object in global memory :\n");
	for (int i = 0; i < 3; i++)
		performanceOnGPU_PoissonGlobalLocal(context, device, queue, "performancePoissonGlobalLocal", n, n1, lambdaArr[0], streamBufferSize, streams, simUseGlobalMem);


	//Free resources
	clrngMrg32k3aDestroyStreams(streams);

	return EXIT_SUCCESS;
}


