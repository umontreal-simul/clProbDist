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
#include "Types.h"
#include "Policies.h"
#include <clRNG/mrg32k3a.h>
#include "../../include/clProbDist/poisson.h"
#include "../../include/clProbDist/gamma.h"


int main()
{

	//double alpha = 9;
	//double lambda = 0.5;
	//int decprec = 15;
	//clprobdistStatus err;

	//clprobdistGamma* gammaDist = clprobdistGammaCreate(alpha, lambda, decprec, &err);
	//clrngMrg32k3aStream* stream = clrngMrg32k3aCreateStreams(NULL, 1, NULL, NULL);

	//double r;

	//for (unsigned i = 0; i < 10; i++){
	//	double u = clrngMrg32k3aRandomU01(stream);
	//	clprobdistGammaInverseCDF(u, gammaDist, &r);

	//	printf("u %f r %f \n", u, r);
	//}


	/*size_t distBufSize;
	clprobdistStatus err;
	double lambda = 520;
	clprobdistPoisson* dist = clprobdistPoissonCreate(lambda, &distBufSize, &err);


	//Generate Poisson Random Numbers
	clrngMrg32k3aStream* stream = clrngMrg32k3aCreateStreams(NULL, 1, NULL, NULL);
	double u;
	int result;
	for (unsigned i = 0; i < 10; i++){
	u = clrngMrg32k3aRandomU01(stream);
	clprobdistPoissonInverseCDFWithObject(dist, u, &result);

	printf("u %20.19f r %d \n", u, result);
	}*/

	int selectedDist = 0, selectedSim = 0;

	printf(" Please select the distribution to use for generating the inventory's demand : \n\n");
	printf(" (1) - Poisson dist. with fixed Lambda.                              \n");
	printf(" (2) - Poisson dist. with random Lambda that follow a Gamma dist. for each day.\n");
	printf(" (3) - Normal dist. truncated at 0\n");
	printf("\n Choice (1 to 3) : ");
	scanf(" %d", &selectedDist); getchar();
	//selectedDist = 3;

	printf("\n\n Please select the simulation to run : \n\n");
	printf(" One Policy : \n");
	printf(" (1) - n runs on CPU, using one stream                                     \n");
	printf(" (2) - n runs with n work items, using two streams and their substreams    \n");
	printf(" (3) - n runs with n work items, using two arrays of n stream each         \n");
	printf(" (4) - n2 runs with n1 work itmes, using 2*n1 streams their substreams     \n");
	printf(" (5) - n2 runs with n1 work items, using 2*n2 distinct streams.            \n");
	printf("\n Several policies  : \n");
	printf(" (6) - p policies in series, with n1 work items and n2 runs per work item\n"); /*option1*/
	printf(" (7) - p policies in parallel, with n1p work items and n2 runs per work item \n"); /*option2*/

	printf("\n Choice (0 to quit) : ");
	scanf(" %d", &selectedSim); getchar();
	//selectedSim = 2;

	int m = 30;  // number of days
	int n =  1 << 10;    // number of runs
	int n1 = 1 << 8;   // number of workitems that will simulate n2 runs 

	double lambda = 50;
	double mu = 100;
	double sigma = 10;

	//Policies values
	int s[] = { 80, 80 };
	int S[] = { 200, 198 };
	int P = 2;  // number of policies

	cl_device_type device_type = CL_DEVICE_TYPE_CPU;
	int platform_index = 0;
	int device_index = 0;
	int ret = 0;

	//Simulate different options :
	while (selectedSim != 0)
	{
		switch (selectedSim)
		{
		case 1: {
			printf("\n===========================================================================\n ONE POLICY:\n=============\n");
			printf("+++++++++     On CPU (basic case): One policy, one stream \n");
			OnePolicyData data = { n, 0, m, s[0], S[0], selectedDist, lambda, mu, sigma, basic };
			ret = call_with_opencl(platform_index, device_type, device_index, &one_Policy, &data, true);

			break;
		}
		case 2: {
			printf("\n===========================================================================\n ONE POLICY:\n=============\n");
			printf("+++++++++     On CPU (case a) : One policy, two streams with their substreams  \n");
			OnePolicyData data = { n, 0, m, s[0], S[0], selectedDist, lambda, mu, sigma, Case_a };
			ret = call_with_opencl(platform_index, device_type, device_index, &one_Policy, &data, true);
			break;
		}
		case 3: {
			printf("\n===========================================================================\n ONE POLICY:\n=============\n");
			printf("+++++++++     On CPU (case b): One policy, two arrays of n streams each  \n");
			OnePolicyData data = { n, 0, m, s[0], S[0], selectedDist, lambda, mu, sigma, Case_b };
			ret = call_with_opencl(platform_index, device_type, device_index, &one_Policy, &data, true);
			break;
		}
		case 4: {
			printf("\n===========================================================================\n ONE POLICY:\n=============\n");
			//no CPU implementation : (case c) 
			OnePolicyData data = { n, n1, m, s[0], S[0], selectedDist, lambda, mu, sigma, Case_c };
			ret = call_with_opencl(platform_index, device_type, device_index, &one_Policy, &data, true);
			break;
		}
		case 5: {
			printf("\n===========================================================================\n ONE POLICY:\n=============\n");
			//no CPU implementation : (case d)
			OnePolicyData data = { n, n1, m, s[0], S[0], selectedDist, lambda, mu, sigma, Case_d };
			ret = call_with_opencl(platform_index, device_type, device_index, &one_Policy, &data, true);
			break;
		}
		case 6: {
			printf("\n===========================================================================\n Several policies Option1:\n=========\n");
			//no CPU implementation : Several policies (option 1)
			SeveralPoliciesData data = { n, n1, m, s, S, P, selectedDist, lambda, mu, sigma, Option1 };
			ret = call_with_opencl(platform_index, device_type, device_index, &several_Policies, &data, true);
			break;
		}
		case 7: {
			printf("\n===========================================================================\n Several policies Option2:\n=========\n");
			//no CPU implementation : Several policies (option 2)
			SeveralPoliciesData data = { n, n1, m, s, S, P, selectedDist, lambda, mu, sigma, Option2 };
			ret = call_with_opencl(platform_index, device_type, device_index, &several_Policies, &data, false);
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
