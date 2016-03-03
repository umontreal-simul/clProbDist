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

#include <clProbDist/poisson.h>
#include <clRNG/mrg31k3p.h>

cl_double simulateOneRun(const clprobdistPoisson* dist, clrngMrg31k3pStream* stream)
{
    cl_double output = 0.0;
    for (int i = 0; i < 100; i++) {
        cl_double u = clrngMrg31k3pRandomU01(stream);
        cl_int poisson_variate = clprobdistPoissonInverseCDFWithObject(dist, u, NULL);
        // modify output using poisson_variate
        output += poisson_variate / (100.0 * clprobdistPoissonMeanWithObject(dist, NULL));
    }
    return output;
}

int main()
{
    // prepare
    int replications = 1024;
    clprobdistPoisson* dist = clprobdistPoissonCreate(50.0, NULL, NULL);
    clrngMrg31k3pStream* stream = clrngMrg31k3pCreateStreams(NULL, 1, NULL, NULL);

    // simulate
    cl_double sum = 0.0;
    for (int i = 0; i < replications; i++)
        sum += simulateOneRun(dist, stream);
    printf("The average output is %.3f\n", sum / replications);

    // clean up
    clprobdistPoissonDestroy(dist);
    clrngMrg31k3pDestroyStreams(stream);

    return 0;
}
