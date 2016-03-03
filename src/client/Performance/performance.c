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
#include "../common.h"
#include "Options.h"
#include <clRNG/mrg32k3a.h>
#include "../../include/clProbDist/poisson.h"
#include "../../include/clProbDist/gamma.h"

#include "SimulateRun.h"


int main()
{
	int selectedSim = 0;

	printf(" Please select the option to execute : \n\n");
	printf(" (1) Option1 :  Performance of PoissonDist with different Lambda using objects and static calls \n");
	printf(" (2) Option2 :  Performance of NormalDist with dist on private memory vs contant memory \n");
	printf(" (3) Option3 :  Performance of PoissontDist where Dist Object is copied on local memory vs global memory \n");

	printf("\n Choice (0 to quit) : ");
	scanf(" %d", &selectedSim); getchar();

	cl_device_type device_type = CL_DEVICE_TYPE_CPU;
	int platform_index = 0;
	int device_index = 0;
	int ret = 0;

	int n = 1 << 25;    // number of runs
	int n1 = 1 << 20;   // number of workitems that will simulate n2 runs 
	double lambdaArr[] = { 2.0, 10.0, 1000.0 };

	//selectedSim = 1;

	//Simulate different options :
	while (selectedSim != 0)
	{
		switch (selectedSim)
		{
		case 1: {
			//Performance of PoissonDist with different Lambda using objects and static calls\n");
			OptionData data = { n, n1, lambdaArr };
			ret = call_with_opencl(platform_index, device_type, device_index, &Option1, &data, true);
			break;
		}

		case 2: {
			//Performance of NormalDist with dist on private memory vs contant memory\n");
			OptionData data = { n, n1 };
			ret = call_with_opencl(platform_index, device_type, device_index, &Option2, &data, true);
			break;
		}
		case 3: {
			//Performance of PoissontDist where Dist Object is copied on local memory vs global memory
			OptionData data = { n, n1, &lambdaArr[1] };
			ret = call_with_opencl(platform_index, device_type, device_index, &Option3, &data, true);
			break;
		}
		default:
			break;
		}

		printf("\n Choice (0 to quit) : ");
		scanf(" %d", &selectedSim);
	}


	return ret;
}
