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

#include <clProbDist/poisson.h>
#include <clRNG/mrg31k3p.h>

int main()
{
	cl_double lambda = 50.0;
	clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, NULL, NULL);
	clrngMrg31k3pStream* stream = clrngMrg31k3pCreateStreams(NULL, 1, NULL, NULL);

	for (int i = 0; i < 30; i++) {

		cl_double u = clrngMrg31k3pRandomU01(stream);

		cl_int with_param = clprobdistPoissonInverseCDF(lambda, u, NULL);

		cl_int with_object = clprobdistPoissonInverseCDFWithObject(dist, u, NULL);

		printf("u=%f, with param/obj=%d/%d   %s\n", u, with_param, with_object,
			with_param == with_object ? "" : "<--");
	}


	clprobdistPoissonDestroy(dist);
	clrngMrg31k3pDestroyStreams(stream);

	return 0;
}