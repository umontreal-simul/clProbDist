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

#include <stdlib.h>
#include "Types.h"
#include "Simulate_Runs.h"

int main()
{
	int n = 1 << 15;    // number of runs
	int n1 = 1 << 10;   // number of workitems that will simulate n2 runs 
	
	double r = 0.05;
	double sigma = 0.5;
	double strike = 100;
	double s0 = 100;
	int s = 12;
	
	// Array zeta[0..s] must contain zeta[0]=0.0, plus the s observation times.
	double* zeta = (double*)malloc((s + 1)*sizeof(double));
	zeta[0] = 0.0;
	for (int j = 1; j <= s; j++)
		zeta[j] = (double)j / (double)s;
	
	//Set Device parameters
	cl_device_type device_type = CL_DEVICE_TYPE_CPU;
	int platform_index = 0;
	int device_index = 0;
	int ret = 0;

	//Start Simulation
	AsianOptionData data = { n, n1, r, sigma, strike, s0, s, zeta };
	ret = call_with_opencl(platform_index, device_type, device_index, &SimulateAsianOption, &data, true);

	return ret;
}
