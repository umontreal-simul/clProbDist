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
 *   David Munger <mungerd@iro.umontreal.ca>        (2015)
 *   Pierre L'Ecuyer <lecuyer@iro.umontreal.ca>     (2015)
 *
 */

#include <clProbDist/poisson.clh>
#include <clRNG/mrg31k3p.clh>

#ifndef NULL
#define NULL ((void*)0)
#endif

double simulateOneRun(__global const clprobdistPoisson* dist, clrngMrg31k3pStream* stream)
{
    double output = 0.0;
    for (int i = 0; i < 100; i++) {
        double u = clrngMrg31k3pRandomU01(stream);
        int poisson_variate = clprobdistPoissonInverseCDFWithObject(dist, u, NULL);
        // modify output using poisson_variate
        output += poisson_variate / (100.0 * clprobdistPoissonMeanWithObject(dist, NULL));
    }
    return output;
}

__kernel void my_kernel(
	__global const clprobdistPoisson* dist,
	__global const clrngMrg31k3pHostStream* host_streams,
	__global double* output)
{
    size_t gid = get_global_id(0);
    clrngMrg31k3pStream stream;
    clrngMrg31k3pCopyOverStreamsFromGlobal(1, &stream, &host_streams[gid]);
    output[gid] = simulateOneRun(dist, &stream);
}
