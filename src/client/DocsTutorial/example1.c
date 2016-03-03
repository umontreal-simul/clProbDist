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

#include <clProbDist/normal.h>
#include <math.h> // for sqrt()
#include <stdio.h>

// #define EFFICIENT

void normalConfidenceInterval(cl_int n, cl_double* observations, cl_double level, cl_double* lower, cl_double* upper)
{
    // compute sample average and variance
    cl_double sum = 0.0;
    cl_double sum_squares = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum         += observations[i];
        sum_squares += observations[i] * observations[i];
    }
    cl_double average  = sum / n;
    cl_double variance = (sum_squares - average * sum) / (n - 1);

#ifdef EFFICIENT
    // compute confidence interval
    cl_double half_width = sqrt(variance / n) *
              clprobdistStdNormalInverseCDF(0.5 * (1.0 + level), NULL);
    *lower = average - half_width;
    *upper = average + half_width;
#else
    // compute normal distribution parameters
    cl_double mu       = average;
    cl_double sigma    = sqrt(variance / n);
    // compute confidence interval
    *lower = clprobdistNormalInverseCDF(mu, sigma, 0.5 * (1.0 - level), NULL);
    *upper = clprobdistNormalInverseCDF(mu, sigma, 0.5 * (1.0 + level), NULL);
#endif
}

int main()
{
    cl_double observations[] = {
            0.39370309,  0.45559733,  0.02482335,  0.84412496,  0.162165,    0.16980586,
            0.16837318,  0.21349251,  0.23347176,  0.67321979,  0.82547116,  0.47291367
    };
    cl_double level = 0.95;
    cl_double lower, upper;
    normalConfidenceInterval(sizeof(observations) / sizeof(observations[0]),
                    observations, level, &lower, &upper);
    printf("%g%% confidence interval: [%.2g, %.2g]\n", level * 100, lower, upper);
    return 0;
}
